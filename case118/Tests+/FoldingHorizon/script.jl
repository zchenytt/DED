import Random, JuMP, Gurobi, Statistics
include("src/Settings.jl");
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Mono.jl");
include("src/Mono2.jl");

Δtˈ1h_ratio, T, S = 15/60, 12, 25;
TimeLimit_per_SP = 4 * 60.

for SEEDi_D = 1:40
    rngD = Random.Xoshiro(SEEDi_D); # for input data generation
    Wscnvec = Vector{Matrix{Float64}}(undef, S);
    zmat = similar(WindGen.zmat); # newer matrix is based on the old ones
    WindGen.reveal!(zmat, 0, Random.Xoshiro(0)); # initialization only once
    WindGen.fill_Wscnvec_once!(Wscnvec, zmat, 1, S); # used to build the initial model
    Ratea_K = [1.04, 1.1, 1.06, 1.0, 1.08];
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

    # This is optimizing mono
    mnt2 = Mono2.buildModel(SEEDi_D, Settings.Env(), S, T, Δtˈ1h_ratio,
        Ratea, F, CaD, EVnode, EVLmax, EVL01, EVEmax, EVEini,
        Bxini, Enode, EEmax, EEini, EPminDiv,
        WPmax,Wscnvec,Wnode, LPmax,LmuTup,Ltype,Lnode,
        GPmax,GPref,Gnode
    );
    Mono2.online_replanning!(mnt2, zmat, Wscnvec, S, T, TimeLimit_per_SP, EEmax, EVLmax, EVEmax, Wnode, WPmax)

    WindGen.reveal!(zmat, 0, Random.Xoshiro(0)); # initialization only once
    WindGen.fill_Wscnvec_once!(Wscnvec, zmat, 1, S); # used to build the initial model
    mnt = Mono.buildModel(SEEDi_D, Settings.Env(), S, T, Δtˈ1h_ratio,
        Ratea, F, CaD, EVnode, EVLmax, EVL01, EVEmax, EVEini,
        Bxini, Enode, EEmax, EEini, EPminDiv,
        WPmax,Wscnvec,Wnode, LPmax,LmuTup,Ltype,Lnode,
        GPmax,GPref,Gnode
    ); Mono.online_replanning!(mnt, zmat, Wscnvec, S, T, TimeLimit_per_SP, EEmax, EVLmax, EVEmax, Wnode, WPmax)

    cost = sum(values(mnt.QY))
    cost2 = sum(values(mnt2.QY))
    println("seed = $SEEDi_D, cost = $cost, cost2 = $cost2")
end



