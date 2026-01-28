module Build
import JuMP, PGLib, PowerModels
import ..Generators, ..DR, ..WindCur, ..Balance, ..Reserve1, ..Reserve2, ..Settings, ..Costs

function _1st_stage_of_mst(m::JuMP.Model, t::Int, T::Int, Other, Reserve, Wind, Demand, Load)
    scene_1, tonly = [1], [t]
    t_to_T, n_to_T = t:T, t+1:T
    # ▶ `δ0` links to the copy variable after scene broadcasting
    δ0 = Generators.add_reserve(m, n_to_T, scene_1, Reserve)
    # `p.s` and `p0.s` links to copy variable after truncate `t=begin`
    p  = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_1, Other)
    p0 = Generators.add(Generators._addc_has_ub, Generators._c2i_addition, m, t_to_T, scene_1, Reserve)
    # `r0` lines to the copy variable after scene broadcasting
    r0 = Reserve1.add_r0(m, t, T, scene_1, Reserve, p0.i) # ◀ [to be copied]
    ϖ = WindCur.add(m, tonly, _ -> scene_1, Wind)
    ζ =      DR.add(m, tonly, _ -> scene_1, Demand)
    Balance.add(m, tonly, scene_1, Wind, Demand, Load, p.s, p0.s, ϖ, ζ)
    Reserve1.add(m, n_to_T, scene_1, Reserve, p0.i, δ0)
    Reserve1.add(m, t, T, scene_1, Reserve, r0, δ0)
    p, p0, δ0, ϖ, ζ, r0
end

inner_model(inn, t, T,#=1:S=# s, Reserve, Wind, Demand, Load) = foreach(wait, [Threads.@spawn(inner_model(inn[s], t, T, s, Reserve, Wind, Demand, Load)) for s=s])
function inner_model(m::JuMP.Model, t::Int, T::Int, s, Reserve, Wind, Demand, Load) # ✅
    scene_s, n_to_T = [s], t+1:T
    # copy variables, to be fixed by 1st-stage trial values
    JuMP.@variable(m,  p[t=n_to_T, s=scene_s])
    JuMP.@variable(m, p0[t=n_to_T, s=scene_s])
    JuMP.@variable(m, r0[t=n_to_T, s=scene_s, z=Reserve.Zone, g=eachindex(Reserve.N[z])])
    JuMP.@variable(m, δ0[t=n_to_T, s=scene_s, z=Reserve.Zone, g=eachindex(Reserve.N[z]), u=(0,1)])
    # 2nd-stage polyhedron
    δ = Generators.add(Generators._addc_no_ub, Generators._c2i_subtract, m, n_to_T, scene_s, Reserve)
    Generators.add_δ_le_δ0_constr(m, n_to_T, scene_s, Reserve, δ.c, δ0)
    ϖ = WindCur.add(m, n_to_T, _ -> scene_s, Wind)
    ζ = DR.add(m, n_to_T, _ -> scene_s, Demand)
    Balance.add(m, n_to_T, scene_s, Wind, Demand, Load, p, p0, ϖ, ζ, δ.s) # register m[:a]
    Reserve2.add(m, t,    s, Reserve, r0, δ.i)
    Reserve2.add(m, t, T, s, Reserve, r0, δ.i)
    JuMP.@expression(m, cost, Costs._2nd_stage(m, n_to_T, s, Reserve, Demand, δ.c, ζ))
    JuMP.@objective(m, Min, cost)
    nothing # copy variables were registered inside inn[s], don't exported separately
end
is_optimality_status(#=inn[s]=# m) = JuMP.has_upper_bound(first(m[:a]))
function set_status_to_feasibility(m)
    a = m[:a]
    JuMP.@objective(m, Min, sum(a))
    foreach(JuMP.delete_upper_bound, a)
end
function set_status_to_optimality(m)
    a = m[:a]
    JuMP.@objective(m, Min, m[:cost])
    foreach(x -> JuMP.set_upper_bound(x, 0.), a)
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
    r0 = Reserve1.add_r0(m, t, T, scene_s, Reserve, p0.i)
    δ = Generators.add(Generators._addc_no_ub, Generators._c2i_subtract, m, n_to_T, scene_s, Reserve)
    # pure 1st-stage constrs
    Balance.add(m, [t], scene_s, Wind, Demand, Load, p.s, p0.s, ϖ, ζ)
    Reserve1.add(m, n_to_T, scene_s, Reserve, p0.i, δ0)
    Reserve1.add(m, t, T, scene_s, Reserve, r0, δ0)
    # Other constrs
    Generators.add_δ_le_δ0_constr(m, n_to_T, scene_s, Reserve, δ.c, δ0)
    Balance.add(m, n_to_T, scene_s, Wind, Demand, Load, p.s, p0.s, ϖ, ζ, δ.s)
    Reserve2.add(m, t,    s, Reserve, r0, δ.i)
    Reserve2.add(m, t, T, s, Reserve, r0, δ.i)
    p, p0, δ0, δ, ϖ, ζ, r0
end
add_pibound_cut(s, m::JuMP.Model, l::ReentrantLock, θ, cˈx, Cn) = @lock(l, JuMP.@constraint(m, cˈx + θ[s] ≥ Cn))
function add_pibound_cut(m::JuMP.Model, l::ReentrantLock, θ, s, cˈx, t, T, Other, Reserve, Wind, Demand, Load)
    _C1 = Settings.Model() # 1st-stage `c` + 2nd-stage `1`
    p, p0, δ0, δ, ϖ, ζ, r0 = Ks_polyhedron(_C1, t, T, s, Other, Reserve, Wind, Demand, Load)
    GC.safepoint()
    Costs.add_C1_obj(_C1, t, T, s, Other, Reserve, Demand, p.c, p0.c, δ0, ζ, δ.c)
    yield()
    JuMP.optimize!(_C1)
    JuMP.termination_status(_C1) === JuMP.OPTIMAL || error()
    add_pibound_cut(s, m, l, θ, cˈx, JuMP.objective_bound(_C1))
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

