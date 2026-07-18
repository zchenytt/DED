import Random, JuMP, Gurobi, Statistics
Random.seed!(hash(2));
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Settings.jl");
# include("src/PreModel.jl") # This is to be manually executed, not a real module
include("src/Models.jl");

S = 5000;
tLeap, T = 3, 12; # 15/5min, number of look-ahead decisions
CaD = Case118.get_Case_Dict();
F = Case118.get_PTDF(CaD);
Ratea = 1.0 * Case118.get_ratea(CaD);

Wnode = Case118.get_windnode();
WPmax = Case118.get_windPmax();
WmuMat = WindGen.mu_06;
Wscnvec = WindGen._S(S, 1; WmuMat=WmuMat);

LmuTup = Case118.l06;
Ltype = Case118.get_load_type();
LPmax = Case118.get_load(CaD);
Lnode = Case118.get_loadnode(CaD);

Ggvec = Case118.get_Ggvec(CaD);
Gnode = Case118.get_Gnode(CaD, Ggvec);
GPmax = Case118.get_GPmax(CaD, Ggvec);
KRmp = 1/8 # for each generator, rampingLim = KRmp * GPmax
GClin = Case118.get_GClin(CaD, Ggvec);
KCres = 5.0 # for each generator, Cres/Cgen = KCres

EClin::Float64 = minimum(GClin)/9; # cost coefficient of (P_chg + P_dischg)
Enode = Case118.get_esnode();
EPmax = Case118.get_esPmax();
EEmax = Case118.get_esEmax();

WCcur = 20 * Statistics.mean(GClin)
LCshd = 5 * WCcur

# run code in PreModel.jl at first to decide these (initial) conditions
EEini = [1., 34.];
Agi = [1, 4, 6, 7, 9, 10, 14, 16, 17, 18 ];
Agp = [1.2625, 4.85, 0.1, 2.23, 1.155, 1.95, 2.257, 2.38875, 6.53, 0.094];

# build the monolithic LP
m, Obj2, p0, en0, en, pe, p, r, ϖ, ζ = Models.model_and_solve(S, tLeap, T, F, Ratea,
    Wnode, WPmax, Wscnvec,
    LmuTup, Ltype, LPmax, Lnode, LCshd,
    Gnode, GPmax, GClin,
    Enode, EPmax, EEmax, EEini, EClin,
    KCres, KRmp, WCcur,
    Agi, Agp
);

JuMP.optimize!(m) # solve as a monolithic LP

for s=1:S # see if we need additional MIP model for ES units
    all(x -> abs(x) < 1e-6, JuMP.value.(pe[:, :, 1, s]) .* JuMP.value.(pe[:, :, 0, s])) || @error("s = $s")
end
