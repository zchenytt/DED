import Random, JuMP, Gurobi, Statistics
include("src/Settings.jl");
include("src/Case118.jl");
include("src/WindGen.jl");
# include("src/Mono.jl");
include("src/ADMM.jl");
# include("src/Models.jl");
# include("src/Ben.jl");

# TODO This myBenUpper(primal heuristic) -> Mono(primal+dual) -> ADMM(dual, post analysis) is practicable.

function result_analysis(mst, sub)
    lb = Ben.run_a_round(mst, sub; mstVType='B', addcut=false);
    ub = Ben.get_ub(mst.Xl, sub, mst.ch);
    agap = ub-lb
    rgap = agap / ub
    printstyled("agap = $agap, rgap = $rgap\n"; color = :magenta)
end

const SEEDi, S_to_W_Ratio = 0, 1
const Δtˈ1h_ratio, T, S = 15/60, 12, (S_to_W_Ratio)ADMM.W_bg
const Ratea_K = [1.04, 1.1, 1.06, 1.0, 1.08]
const rngW = Random.Xoshiro(1+SEEDi) # for wind generation
const rngD = Random.Xoshiro(2+SEEDi) # for input data generation
# when generating test cases, still need to have some additional data, then use Random.Xoshiro(scene_index)

CaD = Case118.get_Case_Dict();
F = Case118.get_PTDF(CaD);
Ratea = rand(rngD, Ratea_K) * Case118.MyRateA; # tune the coefficient
Ggvec = Case118.get_Ggvec(CaD);
Gnode = Case118.get_Gnode(CaD, Ggvec);
GPmax = Case118.get_myGPmax(CaD, Ggvec);
GPref = rand(rngD, .2:1e-3:.8, 1+T, length(GPmax));
LmuTup = Case118.l16; # this furnish values along the time axis
Ltype = Case118.get_load_type(rngD); # for each geo location

LPmax = Case118.get_load(CaD); # for each geo location
Lnode = Case118.get_loadnode(CaD); # for each geo location
Wnode = Case118.get_windnode();
WPmax = Case118.get_windPmax();
Wscnvec = WindGen._S(rngW, S, 1; T=T);
EVnode = Case118.EVnode;
EVLmax = Case118.EVLmax; # This is defined in terms of the network side (meaning that the energy demand at the vehicle side is smaller)
EVL01 = Case118.EVL01;
EVEmax = Case118.EVEmax;
EVEini = Case118.EVEini;
Bxini = Case118.Bxini;
Enode = Case118.get_esnode();
EEmax = Case118.get_esEmax();
EEini = Case118.EEini;
EPminDiv = [(rand(rngD, 7.5:.01:8.5), rand(rngD, 7.5:.01:8.5)) for _ = EEmax]; # the Divisor (discharge, charge)
WmuMat = nothing

adev = Settings.Env(S); admv = similar(adev, JuMP.Model); Threads.@threads for s=1:S
    admv[s] = ADMM.Model(
        adev[s], s, S, T, Δtˈ1h_ratio,
        Ratea, F, CaD,
        EVnode, EVLmax, EVL01, EVEmax, EVEini,
        Bxini, Enode, EEmax, EEini, EPminDiv,
        WPmax,Wscnvec,Wnode,
        LPmax,LmuTup,Ltype,Lnode,
        GPmax,GPref,Gnode
    )
    JuMP.set_attribute(admv[s], "TimeLimit", 48.)
end
ADMM.solve_all_and_check(admv)
# initialize 3 parties
λi, xi_x = [fill(0., length(admv[s][:xi])) for s=1:S], [JuMP.value.(admv[s][:xi]) for s=1:S];
push!(xi_x, sum(xi_x)/S); # xi_x[end] = x_com

res_old = Inf
for k = 1:15
    ADMM.update_dual!(λi, xi_x)
    ADMM.upd_xi_objs!(admv, λi, xi_x)
    ADMM.solve_all_and_check(admv)
    ADMM.upd_all_xi!(xi_x, admv)
    ADMM.upd_x!(xi_x, λi)
    res = ADMM.get_residual(xi_x)
    res > res_old && break
    res_old = min(res_old, res)
    println("k = $k, res = $res, λi = $λi")
end
maximum(abs, sum(λi))
ADMM.post_adjust_λi!(λi) # make sure sum(λi) .== 0
ADMM.post_upd_xi_objs!(admv, λi, xi_x)
ADMM.solve_all_and_check(admv)
[JuMP.MOI.get(m, JuMP.MOI.RelativeGap()) for m=admv]
[JuMP.objective_bound(m) for m=admv]
lb_ADMM = sum(JuMP.objective_value, admv) # -116013.86237403593
