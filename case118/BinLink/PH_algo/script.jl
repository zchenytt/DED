import Random, JuMP, Gurobi, Statistics, LinearAlgebra
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Settings.jl");
include("src/Models.jl");
include("src/Ben.jl");

const SEEDi = 2
const ρ = 9000.
const Ratea_K = 1.
const T, S = 12, 127 * 48
const Δtˈ1h_ratio = 1/4
const ch = Channel{Int}()

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
envs = Settings.Env(S); 
sub = similar(envs, Models.SubMIPTy);
Ben.para_build(
    ch, envs, sub, S, T, Δtˈ1h_ratio,
    F, CaD, Ratea,
    Ggvec, Gnode, GPmax, GPref,
    LmuTup, Ltype, LPmax, Lnode,
    Wnode, WPmax, WmuMat, Wscnvec,
    EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
    EEmax, EEini, EPminDiv, Enode
)
const X = similar(sub[1].Xl); # x_common
k = 0
Threads.@threads for n = sub
    JuMP.set_attribute(n.m, "Method", 1)
    Settings.opt_ass_opt(n)
end
tBotΣ = 0.; t0 = time_ns(); while true
    X .= tbot = 0.; for n = sub
        Gurobi.GRBgetdblattrarray(n.o, "X", 0, n.N, n.Xl) # Xs update
        tbot = max(tbot, Settings.getmodeldblattr(n, "Runtime"))
        X .+= n.Xl
    end; tBotΣ += tbot; X ./= S; # x_common update
    for n = sub
        (; Pi, Xl) = n
        Xl .-= X
        @. Pi = Pi + ρ * Xl # multiplier update
    end
    rsd = maximum(LinearAlgebra.norm(n.Xl, Inf) for n=sub)
    t_elps = 1e-9(time_ns()-t0)
    println("k = $k, rsd = $rsd, t_elps=$t_elps, tBotΣ=$tBotΣ")
    rsd > 1e-6 || break
    tBotΣ > 9e2 && break
    k += 1
    Threads.@threads for n = sub
        (; m, Qy, Pi, xs) = n
        JuMP.@objective(m, Min, Qy + Pi'xs + 0.5ρ * sum((i-j)^2 for (i,j)=zip(xs,X)))
        Settings.opt_ass_opt(n)
    end
end
