import Random, JuMP, Gurobi, Statistics

# src file inside the LowLevel_2SLP is to prove that
# 2-stage stochastic LP can be solved very efficiently
# (using single-theta classic Benders decomposition) even if S is large.

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

function fix_bin_vars(n)
    (; m, o, Pi, N, Bv) = n
    JuMP.set_attribute(m, "SolutionLimit", 1)
    Settings.opt_and_ter(n)
    Bx = similar(Bv, Float64); Gurobi.GRBgetdblattrarray(o, "X", N, length(Bv), Bx)
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv), Bv)
    Gurobi.GRBsetdblattrarray(o, "LB", N, length(Bv), Bx)
    Gurobi.GRBsetdblattrarray(o, "UB", N, length(Bv), Bx)
    JuMP.set_attribute(m, "SolutionLimit", 2000000000)
    Gurobi.GRBreset(o, 0)
    Settings.opt_ass_opt(n)
end
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
    mst = Models.single_θ_mst!(Settings.Env(), S, Δtˈ1h_ratio, EVLmax, EVL01, EVEini, EEmax, EEini, EPminDiv, EVEmax)

    Ben.para_build(
        mst.ch, envs, sub, S, T, Δtˈ1h_ratio,
        F, CaD, Ratea,
        Ggvec, Gnode, GPmax, GPref,
        LmuTup, Ltype, LPmax, Lnode,
        Wnode, WPmax, WmuMat, Wscnvec,
        EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
        EEmax, EEini, EPminDiv, Enode
    )
    Threads.@threads(for n=sub fix_bin_vars(n) end)
    for n=sub, x=n.m[:seS]
        JuMP.set_upper_bound(x, 1e100)
    end
    mst, sub
end
function result_analysis(mst, sub)
    lb = Ben.run_a_round(mst, sub; mstVType='B', addcut=false);
    ub = Ben.get_ub(mst.Xl, sub, mst.ch);
    agap = ub-lb
    rgap = agap / ub
    printstyled("agap = $agap, rgap = $rgap\n"; color = :magenta)
end
single_theta_LP_Ben_train(mst, sub) = for k = 0:999999999
    Settings.opt_ass_opt(mst)
    lb = Settings.getmodeldblattr(mst, "ObjBound")
    Gurobi.GRBgetdblattrarray(mst.o, "X", 1, mst.N, mst.Xl)
    Threads.@threads for n=sub
        Gurobi.GRBsetdblattrarray(n.o, "LB", 0, mst.N, mst.Xl)
        Gurobi.GRBsetdblattrarray(n.o, "UB", 0, mst.N, mst.Xl)
        Settings.opt_ass_opt(n)
        Gurobi.GRBgetdblattrarray(n.o, "RC", 0, mst.N, n.Pi)
        n.Cd[end-1] = obj = Settings.getmodeldblattr(n, "ObjVal")
        n.Cd[1:end-2] .= n.Pi
        n.Cd[end] = Cn = obj - mst.Xl'n.Pi
    end
    mst.Cd .= 0;
    for n=sub mst.Cd .+= n.Cd end
    mst.Cd ./= mst.S;
    ub = mst.Cd[end-1]
    rgap = abs((ub-lb)/ub) 
    println("k = $k, lb = $lb, ub = $ub, rgap = $rgap")
    lb + .8 > ub && break
    mst.Cd[mst.N+1] = -1;
    Gurobi.GRBaddconstr(mst.o, mst.N+1, mst.Ci, mst.Cd, Cchar('<'), -mst.Cd[end], C_NULL)
    k == 0 && Settings.setxdblattrelement(mst, 0, "Obj", 1.)
end

function one_test(SEEDi, a...)
    println("SEEDi = $SEEDi")
    mst, sub = get_mst_sub(a..., SEEDi)
    @time single_theta_LP_Ben_train(mst, sub)
end

test_runtime(a...) = for SEEDi = 0:20
    one_test(SEEDi, a...)
end

const S_to_W_Ratio = 48
const Δtˈ1h_ratio, T, S = 15/60, 12, (S_to_W_Ratio)Ben.W_bg
const Ratea_K = 1.0
test_runtime(Δtˈ1h_ratio, T, S, Ratea_K)

