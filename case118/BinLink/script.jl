import Random, JuMP, Gurobi, Statistics
# include("src/Case118.jl");
# include("src/WindGen.jl");
# include("src/Settings.jl");
# include("src/Models.jl");
# include("src/Ben.jl");

include("Case118.jl");
include("WindGen.jl");
include("Settings.jl");
include("Models.jl");
include("Ben.jl");

function get_mst_sub(Δtˈ1h_ratio, T, S, Ratea_K, SEEDi)
    Random.seed!(hash(SEEDi))
    CaD = Case118.get_Case_Dict();
    F = Case118.get_PTDF(CaD);
    Ratea = Ratea_K * Case118.MyRateA; # tune the coefficient

    Ggvec = Case118.get_Ggvec(CaD);
    Gnode = Case118.get_Gnode(CaD, Ggvec);
    GPmax = Case118.get_myGPmax(CaD, Ggvec);
    GPref = rand(.2:1e-3:.8, 1+T, length(GPmax));

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
    mst = Models.multi_θ_mst!(Settings.Env(), S, EVLmax, EEmax);
    Ben.para_build(
        mst.ch, envs, sub, S, T, Δtˈ1h_ratio,
        F, CaD, Ratea,
        Ggvec, Gnode, GPmax, GPref,
        LmuTup, Ltype, LPmax, Lnode,
        Wnode, WPmax, WmuMat, Wscnvec,
        EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
        EEmax, EEini, EPminDiv, Enode
    )
    mst, sub
end
function result_analysis(mst, sub)
    lb = Ben.run_a_round(mst, sub; mstVType='B', addcut=false);
    ub = Ben.get_ub(mst.Xl, sub, mst.ch);
    agap = ub-lb
    rgap = agap / ub
    printstyled("agap = $agap, rgap = $rgap\n"; color = :magenta)
end
function test(Δtˈ1h_ratio, T, S, Ratea_K, SEEDi)
    printstyled("    test> Ratea_K=$Ratea_K, SEEDi=$SEEDi ¶\n"; color = 24)
    mst, sub = @time "Data and JuMP Models" get_mst_sub(Δtˈ1h_ratio, T, S, Ratea_K, SEEDi);
    @time "furnish_ini_lb" Ben.run_a_round(mst, sub; mstVType='C', addcut=true)
    @time "root_train" Ben.root_train(mst, sub, 12 * 60., 135 * 60.);
    @time "leaf_train" Ben.leaf_train(mst, sub)
    @time "result_analysis" result_analysis(mst, sub)
end

const S_to_W_Ratio = 48
const Δtˈ1h_ratio, T, S = 15/60, 12, (S_to_W_Ratio)Ben.W_bg
const Ratea_K = [1.04, 1.1, 1.06, 1.0, 1.08]

for i = 0:1000
    test(Δtˈ1h_ratio, T, S, Ratea_K[1], i)
    circshift!(Ratea_K,-1)
end
