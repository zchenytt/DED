import Random; Random.seed!(hash(5));
import .Threads.@spawn
const t, T, S = 1, 24, 3
foreach(include, ("Settings.jl", "LoadGen.jl", "WindGen.jl"));
const Other, Reserve, Demand, Load = LoadGen.get_case2383(T);
const Wind = WindGen.get_case2383(t, T, S);
foreach(include, ("WindCur.jl", "DR.jl"));
foreach(include, ("Generators.jl", "Reserve1.jl", "Reserve2.jl", "Balance.jl"));
foreach(include, ("Costs.jl","Build.jl","Get.jl"));

# Status on 2026/01/28:
# it appears that the Benders Decomposition can converges (close the Gap)
# TODO add async programming

import JuMP, Gurobi
const inn = Settings.Model(S);
Build.inner_model(inn, t, T, 1:S, Reserve, Wind, Demand, Load)
const mst, AddCutLock, DeRefLock = Settings.Model(), ReentrantLock(), ReentrantLock();
const p, p0, δ0, ϖ, ζ, r0 = Build._1st_stage_of_mst(mst, t, T, Other, Reserve, Wind, Demand, Load);
const x_link = (p=view(p.s,t+1:T,[1]), p0=view(p0.s,t+1:T,[1]), r0=r0, δ0=δ0);
const cˈx = Costs._1st_stage(mst, t, T, #=scene=# 1, Other, Reserve, Demand, p.c, p0.c, δ0, ζ);
const θ = Costs.add_θ_set_obj(mst, S, cˈx);
foreach(wait, [@spawn(Build.add_pibound_cut(mst, AddCutLock, θ, s, cˈx, t, T, Other, Reserve, Wind, Demand, Load)) for s=1:S])

print_vio(s, str, vio) = @ccall(printf("s=%d-%s vio = %e\n"::Cstring; s::Cint, str::Cstring, vio::Cdouble)::Cint)
function scene_subprocedure(m, s, AddCutLock, DeRefLock, R, x_link, θs::JuMP.VariableRef, mst)
    n = @lock(DeRefLock, R[])
    θ, X = n.θ[s], n.x # trial solution about this_scene
    properties = (:p, :p0, :r0, :δ0)
    Get.fix(m, X, properties)
    JuMP.optimize!(m)
    st = JuMP.termination_status(m)
    if st === JuMP.OPTIMAL
        str, vio, exp = "Opt.", -θ, -θs
    elseif st === JuMP.INFEASIBLE
        Build.set_status_to_feasibility(m)
        JuMP.optimize!(m)
        JuMP.termination_status(m) === JuMP.OPTIMAL || error(4392857)
        str, vio, exp = "Fea.", 0, JuMP.AffExpr(0)
    end
    obj = JuMP.objective_value(m)
    vio += obj
    #=Debug=# print_vio(s, str, vio)
    if vio > 1e-6
        NegCn = Get.innerstats!(#=mutated=# exp, -obj, m, x_link, X, properties)
        @lock(AddCutLock, JuMP.@constraint(mst, exp ≤ NegCn))
    end
    Build.is_optimality_status(m) || Build.set_status_to_optimality(m)
    nothing
end
JuMP.optimize!(mst)
JuMP.termination_status(mst) === JuMP.OPTIMAL || error()
const R = Ref(Get.trial_solution_nt(mst, θ, t+1:T, p.s, p0.s, r0, δ0));

function main()
    for k = 1:300
        for s = 1:S
            scene_subprocedure(inn[s], s, AddCutLock, DeRefLock, R, x_link, θ[s], mst)
        end
        JuMP.optimize!(mst)
        JuMP.termination_status(mst) === JuMP.OPTIMAL || error()
        n = Get.trial_solution_nt(mst, θ, t+1:T, p.s, p0.s, r0, δ0) # This must be allocated anew!
        lb = n.lb
        println("k=$k, lb=$(lb)")
        @lock(DeRefLock, setfield!(R, :x, n));
    end
end
main()
