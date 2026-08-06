import Random, JuMP, Gurobi, Statistics
const SEEDi = 7;
Random.seed!(hash(SEEDi));
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Settings.jl");
include("src/Models.jl");
include("src/Ben.jl");

const Δtˈ1h_ratio, T, S = 15/60, 12, 100;
CaD = Case118.get_Case_Dict();
F = Case118.get_PTDF(CaD);
Ratea = 1.0 * Case118.MyRateA; # tune the coefficient

# generators (e.g. thermal/hydro plants that are re-dispatchable)
Ggvec = Case118.get_Ggvec(CaD);
Gnode = Case118.get_Gnode(CaD, Ggvec);
GPmax = Case118.get_myGPmax(CaD, Ggvec);
GPref = rand(.2:1e-3:.8, 1+T, length(GPmax));

# regular Load
LmuTup = Case118.l16; # this furnish values along the time axis
Ltype = Case118.get_load_type(); # for each geo location
LPmax = Case118.get_load(CaD); # for each geo location
Lnode = Case118.get_loadnode(CaD); # for each geo location

Wnode = Case118.get_windnode();
WPmax = Case118.get_windPmax();
WmuMat = WindGen.mu_06;
Wscnvec = WindGen._S(S, 1; T=T, WmuMat=WmuMat);

EVnode = Case118.EVnode;
EVLmax = Case118.EVLmax; # This is defined in terms of the network side (meaning that the energy demand at the vehicle side is smaller)
EVL01 = Case118.EVL01;
EVEmax = Case118.EVEmax;
EVEini = Case118.EVEini;
Bxini = Case118.Bxini;

Enode = Case118.get_esnode();
EEmax = Case118.get_esEmax();
EEini = Case118.EEini;
EPminDiv = [(rand(7.5:.01:8.5), rand(7.5:.01:8.5)) for _ = EEmax]; # the Divisor (discharge, charge)

envs = Settings.Env(S); sub = similar(envs, Models.SubMIPTy);
Threads.@threads for s=eachindex(envs)
    Models.MIP_sub!(
        sub, s, S, envs[s], T, Δtˈ1h_ratio,
        F, CaD, Ratea,
        Ggvec, Gnode, GPmax, GPref,
        LmuTup, Ltype, LPmax, Lnode,
        Wnode, WPmax, WmuMat, Wscnvec,
        EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
        EEmax, EEini, EPminDiv, Enode
    );
end;
mst = Models.multi_θ_mst!(Settings.Env(), S, EVLmax, EEmax);
Ben.run_a_round(mst, sub; mstVType='C', addcut=true)
ta = Threads.@spawn(:interactive, Ben.root_train(mst, sub, 12 * 60., 7 * 60.));

# ----- early stop -----
Threads.atomic_and!(mst.proceed, false);
wait(ta); Threads.atomic_or!(mst.proceed, true);
# ----------------------

Ben.leaf_train(mst, sub)
lb = Ben.run_a_round(mst, sub; mstVType='B', addcut=false)
ub = Ben.get_ub(mst.Xl, sub);

# x_vec = Vector{Float64}[];
# ub_vec = fill(NaN, 2^6);
# for x in Iterators.product(fill((0., 1.), 6)...)
#     push!(x_vec, collect(x))
# end
# for (i, Xl)=enumerate(x_vec)
#     all(==(1), Xl) || continue
#     Threads.@threads for s=eachindex(sub)
#         Ben.get_ub(s, Xl, sub)
#     end
#     ubv = [n.Cd[end] for n=sub]
#     ub_vec[i] = ub_i = Statistics.mean(ubv)
#     errsum = sum(n.Cd[end-1] for n=sub)
#     println("trial_i = $i, ubv=$ubv, ub_i = $ub_i, errsum = $errsum")
# end
