module Build
import JuMP, Gurobi #, PGLib, PowerModels
import ..Generators, ..DR, ..WindCur, ..Balance, ..Reserve1, ..Reserve2, ..Settings, ..Costs

_3(i) = Vector{Float64}(undef, i)
function get_nt(m)
    xlen, S, o = m[:xlinklen], length(m[:θ]), m.moi_backend
    xlink, θ = _3(xlen), _3(S)
    Gurobi.GRBgetdblattrarray(o, "X", m[:xlinkstart], xlen,    xlink)
    Gurobi.GRBgetdblattrarray(o, "X", m[:θstart],     Cint(S), θ)
    m[:k] = k = m[:k] + 1
    (
        k = k,
        lb = Settings.getmodeldblattr(m, "ObjBound"),
        θ = θ,
        cˈx = Settings.getxdblattrelement(m, m[:cˈx], "X"),
        xlink = xlink
    )
end

_4(m, x, v) = for x=x Settings.setxdblattrelement(m, x, "UB", v) end
function chgabnd(n, ub) # ub=1: Feas. mode; ub=0: minC. mode
    u = Cdouble(ub)
    _4(n, n[:a], 1e20 * u)
    Settings.setxdblattrelement(n, n[:asum], "Obj", u)
    Settings.setxdblattrelement(n, n[:csum], "Obj", 1-u)
    n[:ismC] = !Bool(ub)
    nothing
end

function _mst(m::JuMP.Model, t::Int, T::Int, S, Other, Reserve, Wind, Demand, Load) # ✅
    o = JuMP.backend(m)
    scene_1, tonly, t_to_T, n_to_T = [1], [t], t:T, t+1:T
    # `p.s` and `p0.s` links to copy variable after truncate `t=begin`
    m[:p] = p = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_1, Other)
    m[:p0] = p0 = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_1, Reserve)
    ### Note: These three variables has contiguous ordered Gurobi indices!########
    pp = Generators.add_pp(m, t_to_T, scene_1, p.s, p0.s)
    # `r0` links to the copy variable after scene broadcasting
    m[:r0] = r0 = Reserve1.add_r0(m, t, T, scene_1, Reserve, p0.i) # ◀ [to be copied]
    # ▶ `δ0` links to the copy variable after scene broadcasting
    m[:δ0] = δ0 = Generators.add_reserve(m, n_to_T, scene_1, Reserve) # [t, s, z, g, u]
    ##############################################################################
    m[:ϖ] = ϖ = WindCur.add(m, tonly, _ -> scene_1, Wind)
    m[:ζ] = ζ = DR.add(m, tonly, _ -> scene_1, Demand)
    # `pp` here is full 
    Balance.add1(m, tonly, scene_1, Wind, Demand, Load, pp, ϖ, ζ)
    # `pp` here has the same length to the 2nd-stage subproblem
    m[:pp] = pp2 = pp[n_to_T, scene_1]
    m[:xlinkstart] = Settings._gcc(o, first(pp2))
    m[:xlinklen] = Cint(length(pp2) + length(r0) + length(δ0))
    Reserve1.add6(m, n_to_T, scene_1, Reserve, p0.i, δ0)
    Reserve1.add7(m, t, T, scene_1, Reserve, r0, δ0)
    JuMP.@variable(m, cˈx)
    JuMP.@constraint(m, cˈx == Costs._1st_stage(m, t, T, scene_1, Other, Reserve, Demand, p.c, p0.c, δ0, ζ))
    JuMP.@variable(m, θ[1:S])
    m[:θstart] = Settings._gcc(o, first(θ))
    JuMP.@variable(m, θsum)
    JuMP.@constraint(m, θsum == sum(θ))
    JuMP.@objective(m, Min, cˈx + 1/S * θsum)
    m[:r] = Ref{Cdouble}()
    m[:k] = -1 # label of nt
    m[:Ncuts] = 0
    nothing
end
function inner_model(tks, inn, a...)
    for (s, n)=enumerate(inn)
        tks[s] = Threads.@spawn(inner_model(n, s, a...))
    end
    foreach(wait, tks)
end
function fix(n, xtrial)
    len = Cint(length(xtrial))
    o, st = JuMP.backend(n), n[:xlinkstart]
    Gurobi.GRBsetdblattrarray(o, "LB", st, len, xtrial)
    Gurobi.GRBsetdblattrarray(o, "UB", st, len, xtrial)
end
function inner_model(n::JuMP.Model, s, t::Int, T::Int, Reserve, Wind, Demand, Load, mst) # ✅
    o = JuMP.backend(n)
    scene_s, n_to_T = [s], t+1:T
    # copy variables, to be fixed by 1st-stage trial values
    JuMP.@variable(n, pp[t=n_to_T, s=scene_s]) # (p + p0)[t, s]
    JuMP.@variable(n, r0[t=n_to_T, s=scene_s, z=Reserve.Zone, g=eachindex(Reserve.N[z])])
    JuMP.@variable(n, δ0[t=n_to_T, s=scene_s, z=Reserve.Zone, g=eachindex(Reserve.N[z]), u=(0,1)])
    n[:xlinkstart] = Settings._gcc(o, first(pp))
    # 2nd-stage polyhedron
    n[:δ] = δ = Generators.add(Generators._addc_no_ub, Generators._c2i_subtract, n, n_to_T, scene_s, Reserve)
    Generators.add_δ_le_δ0_constr(n, n_to_T, scene_s, Reserve, δ.c, δ0)
    n[:ϖ] = ϖ = WindCur.add(n, n_to_T, _ -> scene_s, Wind)
    n[:ζ] = ζ = DR.add(n, n_to_T, _ -> scene_s, Demand)
    Balance.add2(n, n_to_T, scene_s, Wind, Demand, Load, pp, ϖ, ζ, δ.s) # register n[:ismC], n[:asum], n[:a]
    Reserve2.add(n, t,    s, Reserve, r0, δ.i)
    Reserve2.add(n, t, T, s, Reserve, r0, δ.i)
    JuMP.@variable(n, csum)
    JuMP.@constraint(n, csum == Costs._2nd_stage(n, n_to_T, s, Reserve, Demand, δ.c, ζ))
    JuMP.set_objective_sense(n, JuMP.MIN_SENSE)
    Settings.setxdblattrelement(n, csum, "Obj", 1.)
    n[:r] = Ref{Cdouble}()
    n[:k] = -1 # if nt.k !== n[:k], then go to fix and n[:k] = nt.k
    _7(mst, s, n)
    nothing
end
function _7(mst, s, n) # help to construct a Benders' cut
    is = Settings._gcc(JuMP.backend(mst), mst[:θ][s])
    ci = Cint[range(mst[:xlinkstart]; length = mst[:xlinklen]); is]
    n[:cd] = similar(ci, Cdouble)
    n[:ci] = ci # read-only
    # apart from these, write `cn` to Gurobi directly
    # remember to write `cd[end]` to `0.0` or `-1.0`
    # remember to notice the sign `+/-` of the coefficients and `cn`
end

function Ks_polyhedron(m::JuMP.Model, t::Int, T::Int, s, Other, Reserve, Wind, Demand, Load)
    scene_s, t_to_T, n_to_T = [s], t:T, t+1:T
    # 1st-stage decision: p, p0, δ0, ϖ(t), ζ(t)
    # 2nd-stage decision: δ(all), ϖ(τ:T), ζ(τ:T)
    ϖ = WindCur.add(m, t_to_T, _ -> scene_s, Wind)
    ζ =      DR.add(m, t_to_T, _ -> scene_s, Demand)
    δ0 = Generators.add_reserve(m, n_to_T, scene_s, Reserve)
    p = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_s, Other)
    p0 = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_s, Reserve)
    pp = Generators.add_pp(m, t_to_T, scene_s, p.s, p0.s)
    r0 = Reserve1.add_r0(m, t, T, scene_s, Reserve, p0.i)
    δ = Generators.add(Generators._addc_no_ub, Generators._c2i_subtract, m, n_to_T, scene_s, Reserve)
    # pure 1st-stage constrs
    Balance.add1(m, [t], scene_s, Wind, Demand, Load, pp, ϖ, ζ)
    Reserve1.add6(m, n_to_T, scene_s, Reserve, p0.i, δ0)
    Reserve1.add7(m, t, T, scene_s, Reserve, r0, δ0)
    # Other constrs
    Generators.add_δ_le_δ0_constr(m, n_to_T, scene_s, Reserve, δ.c, δ0)
    Balance.add2(m, n_to_T, scene_s, Wind, Demand, Load, pp, ϖ, ζ, δ.s)
    Reserve2.add(m, t,    s, Reserve, r0, δ.i)
    Reserve2.add(m, t, T, s, Reserve, r0, δ.i)
    p, p0, δ0, δ, ϖ, ζ, r0
end

function add_pibound_cut(tks, a...)
    for s = eachindex(tks)
        tks[s] = Threads.@spawn(add_pibound_cut(s, a...))
    end
    foreach(wait, tks)
end
_9(m, cx, th, Cn) = JuMP.@constraint(m, cx + th ≥ Cn)
function add_pibound_cut(s::Int, m::JuMP.Model, l::ReentrantLock, t, T, Other, Reserve, Wind, Demand, Load)
    θ, cˈx = m[:θ], m[:cˈx]
    _C1 = Settings.Model() # 1st-stage `c` + 2nd-stage `1`
    p, p0, δ0, δ, ϖ, ζ, r0 = Ks_polyhedron(_C1, t, T, s, Other, Reserve, Wind, Demand, Load)
    Costs.add_C1_obj(_C1, t, T, s, Other, Reserve, Demand, p.c, p0.c, δ0, ζ, δ.c)
    JuMP.optimize!(_C1)
    JuMP.termination_status(_C1) === JuMP.OPTIMAL || error()
    th, Cn = θ[s], JuMP.objective_bound(_C1)
    @lock(l, _9(m, cˈx, th, Cn))
    nothing
end;


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
