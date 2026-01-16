import Random; Random.seed!(hash(11));
foreach(include, ("LoadGen.jl", "WindGen.jl", "Dispatch.jl"));
const T, S = 24, 10000; # 5000 scenarios already reaches the 1 hour time limit
const Z, F, B, N, FlowLimitNameplate = Dispatch.load_network();
initiate(GRB_ENV) = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV));
Demand, Load, Other, Reserve = LoadGen.get_case2383(T);
Wind = WindGen.get_case2383(S); 

import JuMP, Gurobi
const GRB_ENV = Gurobi.Env();
const model = initiate(GRB_ENV);
ϖ = Dispatch.add_ϖ(model, T, S, Wind);
ζ = Dispatch.add_ζ(model, T, S, Demand);
p = Dispatch.add_power_and_c2i_constrs(model, Other, T);
p0 = Dispatch.add_power_and_c2i_constrs(model, Reserve, T);
δ0 = Dispatch.add_δ0(model, Reserve, T);
δ = Dispatch.add_δ(model, Reserve, T, S);

# pure 1st-stage
Dispatch.add_power_balance_1(model, Other, Reserve, Wind, Load, Demand, p.i, p0.i, ϖ, ζ);
Dispatch.add_1st_stage_reserve_constrs(model, Reserve, T, δ0, p0.i);

# 1st-2nd linking constrs
Dispatch.add_δ_le_δ0_constr(model, Reserve, T, S, δ, δ0);
Dispatch.add_ramp_on_realization_constrs(model, Reserve, T, S, p0.i, δ);
Dispatch.add_power_balance_2(model, T, S, Other, Reserve, Wind, Load, Demand, p.i, p0.i, ϖ, ζ, δ);

build_cost(model, T, S, Other, Reserve, Demand, p, p0, δ, δ0, ζ) = (
    p = Dispatch.build_pCost(model, Other, T, p.c),
    p0 = Dispatch.build_pCost(model, Reserve, T, p0.c),
    δ = Dispatch.build_δ_UpCost(model, Reserve, T, S, δ),
    δ0 = Dispatch.build_δ0Cost(model, Reserve, T, δ0),
    ζ1 = Dispatch.build_ζ_1st_stage_cost(model, Demand, ζ),
    ζ2 = Dispatch.build_ζ_2nd_stage_cost(model, T, S, Demand, ζ)
)
function add_obj(model, T, S, Other, Reserve, Demand, p, p0, δ, δ0, ζ)
    Cost = build_cost(model, T, S, Other, Reserve, Demand, p, p0, δ, δ0, ζ)
    JuMP.@objective(model, Min, sum(i for i = Cost))
    Cost
end;
Cost = add_obj(model, T, S, Other, Reserve, Demand, p, p0, δ, δ0, ζ);

JuMP.set_attribute(model, "Crossover", 0);
JuMP.set_attribute(model, "Method", 2);
# JuMP.set_attribute(model, "PDHGGPU", 1);
# JuMP.set_attribute(model, "TimeLimit", 3600);
JuMP.optimize!(model)

