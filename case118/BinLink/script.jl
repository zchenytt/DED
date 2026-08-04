import Random, JuMP, Gurobi, Statistics
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Settings.jl");
include("src/Models.jl");
include("src/Ben.jl");

# 🟥 This generate subMIPs with distinct level of hardness
const SEEDi = 4; # 5 is very hard MIP
const T, S = 12, 310;
const Δtˈ1h_ratio = 15/60;
Random.seed!(hash(SEEDi));
CaD = Case118.get_Case_Dict();
F = Case118.get_PTDF(CaD);
Ratea = 1.0 * Case118.MyRateA; # tune the coefficient

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

LmuTup = Case118.l16; # this furnish values along the time axis
Ltype = Case118.get_load_type(); # for each geo location
LPmax = Case118.get_load(CaD); # for each geo location
Lnode = Case118.get_loadnode(CaD); # for each geo location

Enode = Case118.get_esnode();
EEmax = Case118.get_esEmax();
EEini = Case118.EEini;
EPminDiv = [(rand(7.5:.01:8.5), rand(7.5:.01:8.5)) for _ = EEmax]; # the Divisor (discharge, charge)

envs = Settings.Env(S); sub = similar(envs, Models.SubMIPTy);
Threads.@threads for s=eachindex(envs)
    Models.MIP_sub!(
        sub, s, S, envs[s], T, Δtˈ1h_ratio,
        F, CaD, Ratea,
        Wnode, WPmax, WmuMat, Wscnvec,
        LmuTup, Ltype, LPmax, Lnode,
        EEmax, EEini, EPminDiv, Enode,
        EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
    );
end;
mst = Models.multi_θ_mst!(Settings.Env(), S, EVLmax, EEmax);
Ben.MIP_warm(mst, sub);
Ben.main(mst, sub, 12 * 60., 3.4 * 60.)
Ben.get_lb_ub(mst, sub)
