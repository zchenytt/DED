# According to experiments results:
# with pure Continuous linking variables, the resulting gap is high (perhaps due to high ub)
# had best follow the setting in BDD2020

# For generating SB cut from MIP subproblems, the TimeLimit should be tuned carefully, based on the general hardness of MIPs

import Random, JuMP, Gurobi, Statistics
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Settings.jl");
include("src/Models.jl");
include("src/Ben.jl");
include("src/MultiTheta.jl");

const T, S = 12, 1000
function _12(SEEDi)
    Random.seed!(hash(SEEDi));
    CaD = Case118.get_Case_Dict();
    F = Case118.get_PTDF(CaD);
    Ratea = 1.0 * Case118.get_ratea(CaD);

    Wnode = Case118.get_windnode();
    WPmax = Case118.get_windPmax();
    WmuMat = WindGen.mu_06;
    Wscnvec = WindGen._S(S, 1; T=T, WmuMat=WmuMat);

    LmuTup = Case118.l16;
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
    EPmin = Case118.get_esPmin(); # 🟧 For MIP ES model only; There're mainly 2 motivations for MIP ES, one is this nonzeroness, the other is simultaneous ch/dis
    EEmax = Case118.get_esEmax();

    WCcur = 20 * Statistics.mean(GClin);
    LCshd = 5 * WCcur;

    # run code in PreModel.jl at first to decide these (initial) conditions
    EEini = [.35, .75] .* EEmax;
    Agi = [1, 4, 6, 7, 9, 10, 14, 16, 17, 18 ];
    Agp = [1.2625, 4.85, 0.1, 2.23, 1.155, 1.95, 2.257, 2.38875, 6.53, 0.094];

    genv = Settings.Env();
    mst = Models.multi_θ_mst!(genv, S, Agi, EEmax);
    sub = Vector{Models.SubMIPTy}(undef, S); Threads.@threads for s=1:S
        Models.MIP_sub!(
            sub, s, S, mst.N, T, F, genv, Agi, Agp, Ratea,
            Lnode, Ltype, LmuTup, LPmax, LCshd,
            WPmax, Wscnvec, Wnode, WCcur,
            EEmax, EPmax, EEini, Enode, EPmin,
            KRmp, GPmax, Gnode, GClin, KCres,
        )
    end
    mst, sub
end;
mst, sub = _12(3);
@time "MultiTheta.MIP_warm (Execute only once)" MultiTheta.MIP_warm(mst, sub)

WC_TIME = 15 * 60.

for i = 1:600
    MultiTheta.train_asyncly_shortly(mst, sub, WC_TIME)
    MultiTheta.eval_at_sync_point(mst, sub, 'B')
end

