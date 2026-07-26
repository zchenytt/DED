import Random, JuMP, Gurobi, Statistics

printstyled("script.jl>\n"; color = 36)

# Use these in REPL
# include("src/Case118.jl");
# include("src/WindGen.jl");
# include("src/Settings.jl");
# include("src/Models.jl");
# include("src/Ben.jl");

# Use these in Bash
include("Case118.jl");
include("WindGen.jl");
include("Settings.jl");
include("Models.jl");
include("Ben.jl");

function mymain(train_seconds)
    Random.seed!(hash(3));
    S = 10000;
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

    genv = Settings.Env()
    mst = Models.mst!(genv, S, Agi, Agp, GPmax, KRmp, EEmax, EEini, EPmax);
    sub = Vector{Models.SubTy}(undef, S); Ben.sub!( # This builds in parallel
        sub, S, mst.N, T, F, tLeap, genv, Agi, Ratea,
        Lnode, Ltype, LmuTup, LPmax, LCshd,
        WPmax, Wscnvec, Wnode, WCcur,
        EEmax, EPmax, EEini, Enode, EClin,
        KRmp, GPmax, Gnode, GClin, KCres
    )
    Settings.opt_ass_opt(mst)
    Ben.load_che(mst)
    @time "warm" Ben.warm(mst, sub, S)
    Ben.set_mst_θ_obj(mst)
    Settings.opt_ass_opt(mst)
    Ben.eval_gap(mst, sub, S)
    Ben.main(mst, sub, S, train_seconds)
    Ben.eval_gap(mst, sub, S)
    free_GiB = Sys.free_memory() / 1024^3
    println("free_GiB = $free_GiB")
end

mymain(5.)
mymain(10.)
mymain(15.)
mymain(20.)
mymain(30.)
mymain(40.)


amd@amd:~/julia_projects/ded/c118$ julia --project=. --threads=255,1 src/script.jl
script.jl>
master built with N = 12
warm: 14.937243 seconds (417.40 k allocations: 21.445 MiB, 313 lock conflicts, 150.13% compilation time)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 2697.1390910090445, ub = 5045.614466088944, rgap = 0.4654488349959675, agap = 2348.4753750798995
free_GiB = 216.38591766357422
master built with N = 12
warm: 12.769560 seconds (69.75 k allocations: 4.304 MiB, 540 lock conflicts)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 2741.4077368507124, ub = 3097.2496501067667, rgap = 0.11488964515462548, agap = 355.84191325605434
free_GiB = 140.08990478515625
master built with N = 12
warm: 12.800128 seconds (69.75 k allocations: 4.304 MiB, 521 lock conflicts)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 3075.6361412387505, ub = 3307.5521102128673, rgap = 0.07011710208828462, agap = 231.91596897411682
free_GiB = 83.9359130859375
master built with N = 12
warm: 12.791393 seconds (69.75 k allocations: 4.304 MiB, 546 lock conflicts)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 3084.3480061722903, ub = 3203.194080798978, rgap = 0.037102364586364324, agap = 118.84607462668782
free_GiB = 46.285091400146484
master built with N = 12
warm: 12.735841 seconds (69.75 k allocations: 4.304 MiB, 566 lock conflicts)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 3084.745961255003, ub = 3085.1752224648653, rgap = 0.00013913673581221597, agap = 0.4292612098624886
free_GiB = 52.18516540527344
master built with N = 12
warm: 12.912030 seconds (69.75 k allocations: 4.304 MiB, 637 lock conflicts)
lb = -3366.6675616787916, ub = 6499.48949516399, rgap = 1.5179895381296937, agap = 9866.157056842781
lb = 3085.0656280720927, ub = 3085.0656967642335, rgap = 2.226602202663732e-8, agap = 6.869214075777563e-5
free_GiB = 28.745594024658203
