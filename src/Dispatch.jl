module Dispatch

import JuMP, PGLib, PowerModels

function load_network()
    PowerModels.silence()
    case = "pglib_opf_case24_ieee_rts.m"
    case = "pglib_opf_case118_ieee.m"
    case = "pglib_opf_case2383wp_k.m"
    Z = PowerModels.make_basic_network(PGLib.pglib(case))
    F = PowerModels.calc_basic_ptdf_matrix(Z)
    B, N = size(F)
    FlowLimitNameplate = map(i -> Z["branch"][string(i)]["rate_a"], 1:B)
    Z, F, B, N, FlowLimitNameplate
end

# Wind curtailment decisions
Wind_PMax_Realized(Wind, s, t, z, w) = Wind.S[s][t,z]Wind.PMax[z][w]
add_ϖ(model, T, S, Wind) = JuMP.@variable(model,
    [t=1:T, s=ifelse(t>1, 1:S, 1:1), z=Wind.Zone, w=eachindex(Wind.N[z])],
    lower_bound = 0, upper_bound = Wind_PMax_Realized(Wind, s, t, z, w)
)

# "Demand" response decision
UB_ζ(PMax, Div, t, z, d, u) = PMax[z][t,d]/Div[u]
function add_ζ(model, T, S, Demand)
    Div = (u = 2, d = 1) # `:d` means cutting Demand
    PMax = Demand.P
    JuMP.@variable(model, [t=1:T, s=ifelse(t>1, 1:S, 1:1), z=Demand.Zone, d=eachindex(Demand.N[z]), u=(:u, :d)],
        lower_bound = 0, upper_bound = UB_ζ(PMax, Div, t, z, d, u)
    )
end

# Generation Power
add_pwl_power(model, GD, T, UB) = JuMP.@variable(model, # GD = Reserve || Other
    [z=GD.Zone, g=eachindex(GD.N[z]), t=1:T, h=(:h, :l)],
    lower_bound = 0, upper_bound = UB[h][z][g]
)
add_power(model, GD, T, ::Symbol) = JuMP.@variable(model,
    [z=GD.Zone, g=eachindex(GD.N[z]), t=1:T],
)
prepare_ub(GD) = (
    h = [[(1-GD.pCost[z].KneeFraction[g])GD.PMax[z][g] for g=eachindex(GD.N[z])] for z=GD.Zone],
    l = [[    GD.pCost[z].KneeFraction[g]GD.PMax[z][g] for g=eachindex(GD.N[z])] for z=GD.Zone]
)
add_pwl_power(model, GD, T) = add_pwl_power(model, GD, T, prepare_ub(GD))
add_power(model, GD, T) = (
    i = add_power(model, GD, T, :identity),
    c = add_pwl_power(model, GD, T)
)
add_c2i_constr(model, GD, T, p) = JuMP.@constraint(model,
    [z=GD.Zone, g=eachindex(GD.N[z]), t=1:T],
    p.i[z,g,t] == p.c[z,g,t,:h] + p.c[z,g,t,:l] # TODO: add a nonzero PMin
)
function add_power_and_c2i_constrs(model, GD, T)
    p = add_power(model, GD, T)
    add_c2i_constr(model, GD, T, p)
    p
end

# Reserve && Re-dispatch
add_δ0(model, Reserve, T) = JuMP.@variable(model,
    [z=Reserve.Zone, r=eachindex(Reserve.N[z]), t=2:T, u=(:u, :d)],
    lower_bound = 0,
    upper_bound = Reserve.Ramp[u][z][r]
)
add_δ(model, Reserve, T, S) = JuMP.@variable(model,
    [s=1:S, z=Reserve.Zone, r=eachindex(Reserve.N[z]), t=2:T, u=(:u, :d)],
    lower_bound = 0
)
add_δ_le_δ0_constr(model, Reserve, T, S, δ, δ0) = JuMP.@constraint(model,
    [s=1:S, z=Reserve.Zone, r=eachindex(Reserve.N[z]), t=2:T, u=(:u, :d)],
    δ[s,z,r,t,u] ≤ δ0[z,r,t,u]
)

# 1st-stage constrs for Reserve Units
function add_δ0Dn_p0_δ0Up_constrs(model, p0, δ0, UB, z, r, t)
    JuMP.@constraint(model, p0 + δ0[z,r,t,:u] ≤ UB)
    JuMP.@constraint(model, δ0[z,r,t,:d] ≤ p0) # TODO: add a nonzero PMin
end
add_δ0Dn_p0_δ0Up_constrs(model, Reserve, T, δ0, p0) = for z=Reserve.Zone, r=eachindex(Reserve.N[z]), t=2:T
    add_δ0Dn_p0_δ0Up_constrs(model, p0[z,r,t], δ0, Reserve.PMax[z][r], z, r, t)
end
function add_1st_stage_ramp_constrs(model, p0τ_m_t, δ0, z, r, τ, RampUp, RampDn)
    JuMP.@constraint(model, p0τ_m_t + δ0[z,r,τ,:u] ≤ RampUp)
    JuMP.@constraint(model, δ0[z,r,τ,:d] - p0τ_m_t ≤ RampDn)
end
add_1st_stage_ramp_constrs(model, Reserve, T, δ0, p0) = for z=Reserve.Zone, r=eachindex(Reserve.N[z]), (t, τ) = zip(1:T, 2:T)
    add_1st_stage_ramp_constrs(model, p0[z,r,τ] - p0[z,r,t], δ0, z, r, τ, Reserve.Ramp.u[z][r], Reserve.Ramp.d[z][r])
end
function add_1st_stage_reserve_constrs(a,b,c,d,e)
    add_δ0Dn_p0_δ0Up_constrs(a,b,c,d,e)
    add_1st_stage_ramp_constrs(a,b,c,d,e)
end

# These are ALL 1st-2nd linking constrs
function _add_ramp_constr(model, pτ_m_pt, RampUp, RampDn)
    JuMP.@constraint(model, pτ_m_pt ≤  RampUp)
    JuMP.@constraint(model, pτ_m_pt ≥ -RampDn)
end
function add_ramp_period_ge2_constr(model, Ramp, p0, δ, z, r, s, t, τ)
    pτ_m_pt = JuMP.@expression(model, p0[z,r,τ]+δ[s,z,r,τ,:u]-δ[s,z,r,τ,:d]-p0[z,r,t]-δ[s,z,r,t,:u]+δ[s,z,r,t,:d]) # ✅
    _add_ramp_constr(model, pτ_m_pt, Ramp.u[z][r], Ramp.d[z][r])
end
add_ramp_period_ge2_constr(model, Reserve, T, S, p0, δ) = for s=1:S, z=Reserve.Zone, r=eachindex(Reserve.N[z]), (t, τ)=zip(2:T, 3:T)
    add_ramp_period_ge2_constr(model, Reserve.Ramp, p0, δ, z, r, s, t, τ)
end
function add_ramp_period_1_to_2_constr(model, Ramp, p0, δ, z, r, s)
    (t, τ) = (1, 2)
    pτ_m_pt = JuMP.@expression(model, p0[z,r,τ]+δ[s,z,r,τ,:u]-δ[s,z,r,τ,:d]-p0[z,r,t]) # ✅
    _add_ramp_constr(model, pτ_m_pt, Ramp.u[z][r], Ramp.d[z][r])
end
add_ramp_period_1_to_2_constr(model, Reserve, S, p0, δ) = for s=1:S, z=Reserve.Zone, r=eachindex(Reserve.N[z])
    add_ramp_period_1_to_2_constr(model, Reserve.Ramp, p0, δ, z, r, s)
end
function add_ramp_on_realization_constrs(model, Reserve, T, S, p0, δ)
    add_ramp_period_ge2_constr(model, Reserve, T, S, p0, δ)
    add_ramp_period_1_to_2_constr(model, Reserve, S, p0, δ)
end

# Power Balance Constrs: 1st and 2nd stage
Load_NamePlate(Load, t, z, l) = Load.P[z][t, l]
Load_NamePlate(Load, t) = sum(Load_NamePlate(Load,t,z,l) for z=Load.Zone for l=eachindex(Load.N[z]))
WindGenSys(Wind, t, s) = sum(Wind_PMax_Realized(Wind, s, t,z,w) for z=Wind.Zone for w=eachindex(Wind.N[z]))
GenSysConst(Wind, Load, Demand, t, s) = WindGenSys(Wind, t, s) - Load_NamePlate(Load, t) - Load_NamePlate(Demand, t)
GenSysExpr(model, GD, p, t) = JuMP.@expression(model, sum(p[z,g,t] for z=GD.Zone for g=eachindex(GD.N[z])))
GenBaseVariable(model, Other, Reserve, Wind, Demand, t, s, p, p0, ϖ, ζ) = JuMP.@expression(
    model, # purely decisions included
    GenSysExpr(model, Other, p, t) + GenSysExpr(model, Reserve, p0, t) +
    sum(ζ[t,s,z,d,:d]-ζ[t,s,z,d,:u] for z=Demand.Zone for d=eachindex(Demand.N[z])) -
    sum(ϖ[t,s,z,w] for z=Wind.Zone for w=eachindex(Wind.N[z]))
)
function add_power_balance_1(model, Other, Reserve, Wind, Load, Demand, p, p0, ϖ, ζ)
    (t, s) = (1, 1)
    JuMP.@constraint(model, 
        GenSysConst(Wind, Load, Demand, t, s) +
        GenBaseVariable(model, Other, Reserve, Wind, Demand, t, s, p, p0, ϖ, ζ)
        == 0
    )
end
add_power_balance_2(model, T, S, Other, Reserve, Wind, Load, Demand, p, p0, ϖ, ζ, δ) = JuMP.@constraint(model,
    [s=1:S, t=2:T], 
    GenRedispatch(model, Reserve, δ, t, s) +
    GenSysConst(Wind, Load, Demand, t, s) +
    GenBaseVariable(model, Other, Reserve, Wind, Demand, t, s, p, p0, ϖ, ζ)
    == 0
)
GenRedispatch(model, Reserve, δ, t, s) = JuMP.@expression(model,
    sum(δ[s,z,r,t,:u]-δ[s,z,r,t,:d] for z=Reserve.Zone for r=eachindex(Reserve.N[z]))
)

# TODO Cost Terms
# TODO
# TODO
prepare_h_cost(GD) = [[GD.pCost[z].Slope[g]GD.pCost[z].Mul[g] for g=eachindex(GD.N[z])] for z=GD.Zone]
build_pCost(model, GD, T, p, H_Cost) = JuMP.@expression(model, sum(
    p[z,g,t,:l]GD.pCost[z].Slope[g] +
    p[z,g,t,:h]H_Cost[z][g]
    for z=GD.Zone for g=eachindex(GD.N[z]) for t=1:T
))
build_pCost(model, GD, T, p) = build_pCost(model, GD, T, p, prepare_h_cost(GD))
build_δ_UpCost(model, GD, T, S, δ, H_Cost) = JuMP.@expression(model, sum(
    δ[s,z,r,t,:u]H_Cost[z][r]
    for s=1:S for z=GD.Zone for r=eachindex(GD.N[z]) for t=2:T
))
build_δ_UpCost(model, GD, T, S, δ) = build_δ_UpCost(model, GD, T, S, δ, prepare_h_cost(GD)/S)
build_δ0Cost(model, Reserve, T, δ0) = JuMP.@expression(model, sum(
    Reserve.δ0Cost[u][z][r]δ0[z,r,t,u]
    for z=Reserve.Zone for r=eachindex(Reserve.N[z]) for t=2:T for u=(:u, :d)
))

build_ζ_2nd_stage_cost(model, T, S, Demand, ζ) = JuMP.@expression(model,
    sum(ζ[t,s,z,d,u]Demand.pCost[u][z][t,d]
    for t=2:T for s=1:S for z=Demand.Zone for d=eachindex(Demand.N[z]) for u=(:u, :d)
)/S)
function build_ζ_1st_stage_cost(model, Demand, ζ)
    t, s = 1, 1
    JuMP.@expression(model,
       sum(ζ[t,s,z,d,u]Demand.pCost[u][z][t,d]
       for z=Demand.Zone for d=eachindex(Demand.N[z]) for u=(:u, :d)
   ))
end



function add_f(model, T, B, S, O, R, W, L, D, F, OtherNode, ReserveNode, WindNode, LoadNode, DemandNode, p, p0, δ, ϖ, ζ, Load)
    f = JuMP.@variable(model, [1:T, 1:B, 1:S])
    JuMP.@constraint(model, [t=2:T, b=1:B, s=1:S], f[t, b, s] ==
        sum(F[b, OtherNode[o]]p[t, o, s] for o = 1:O)
        + sum(F[b, ReserveNode[r]]p2(p0, δ, t, r, s) for r = 1:R)
        + sum(F[b, WindNode[w]]power_after_curtail(ϖ[t, w, s]) for w = 1:W)
        - sum(F[b, LoadNode[l]]Load[t, l] for l = 1:L)
        - sum((power_after_curtail(ζ.d[t, d, s]) + ζ.u[t, d, s])F[b, DemandNode[d]] for d = 1:D)
    )
    # t = 1, deterministic (let s = 1)
    JuMP.@constraint(model, [b=1:B], f[1, b, 1] == 
        sum(F[b, OtherNode[o]]p[1, o, 1] for o = 1:O)
        + sum(F[b, ReserveNode[r]]p0[1, r] for r = 1:R)
        + sum(F[b, WindNode[w]]power_after_curtail(ϖ[1, w, 1]) for w = 1:W)
        - sum(F[b, LoadNode[l]]Load[1, l] for l = 1:L)
        - sum((power_after_curtail(ζ.d[1, d, 1]) + ζ.u[1, d, 1])F[b, DemandNode[d]] for d = 1:D)
    )
    f
end
function add_a(model, f) # This is not involved in the final ED formulation
    T, B, S = size(f)
    a = JuMP.@variable(model, [1:B])
    JuMP.@constraint(model, [t=1:T, b=1:B, s=1:S], a[b] >=  f[t, b, s])
    JuMP.@constraint(model, [t=1:T, b=1:B, s=1:S], a[b] >= -f[t, b, s])
    JuMP.@objective(model, Min, sum(a))
    a
end

# TODO adding branch flow limits
# TODO consider adding a regional reserve lower bound

end

