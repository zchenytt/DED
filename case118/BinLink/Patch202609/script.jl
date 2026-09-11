import Random, JuMP, Gurobi, Statistics
include("src/Settings.jl");
include("src/Case118.jl");
include("src/WindGen.jl");
include("src/Mono.jl");
include("src/Models.jl");
include("src/Ben.jl");

function result_analysis(mst, sub)
    lb = Ben.run_a_round(mst, sub; mstVType='B', addcut=false);
    ub = Ben.get_ub(mst.Xl, sub, mst.ch);
    agap = ub-lb
    rgap = agap / ub
    printstyled("agap = $agap, rgap = $rgap\n"; color = :magenta)
end

const SEEDi, S_to_W_Ratio = 11, 1
const Δtˈ1h_ratio, T, S = 15/60, 12, (S_to_W_Ratio)Ben.W_bg
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

Ben.run_a_round(mst, sub; mstVType='C', addcut=true)
Ben.root_train(mst, sub, 5 * 60., 5 * 60.);
Ben.leaf_train(mst, sub)
result_analysis(mst, sub)

# This is optimizing mono
genv = Settings.Env();
m = Mono.Model(
    genv, S, T, Δtˈ1h_ratio,
    Ratea, F, CaD,
    EVnode, EVLmax, EVL01, EVEmax, EVEini,
    Bxini, Enode, EEmax, EEini, EPminDiv,
    WPmax,Wscnvec,Wnode,
    LPmax,LmuTup,Ltype,Lnode,
    GPmax,GPref,Gnode
)
JuMP.unset_silent(m)
JuMP.set_attribute(m, "TimeLimit", 5 * 60)
JuMP.optimize!(m)
grbObjVal = JuMP.objective_value(m)
ObjBnd = JuMP.objective_bound(m)
grbGap = (grbObjVal - ObjBnd) / grbObjVal
for s=1:S
    l = ifelse(s==1, 0, 1)
    for t=l:T
        for e=eachindex(EEmax), u=(0,1)
            v = JuMP.value(sub[s].m[:bES][e, t, u])
            JuMP.fix(m[:bES][e,t,u,s], round(v); force=true)
        end
        for a=eachindex(EVLmax)
            v = JuMP.value(sub[s].m[:bx][a, t])
            JuMP.fix(m[:bx][a, t, s], round(v); force=true)
        end
    end
end
JuMP.set_silent(m)
JuMP.optimize!(m)
JuMP.termination_status(m) == JuMP.OPTIMAL || error()
localObjVal = JuMP.objective_value(m)
localGap = (localObjVal - ObjBnd) / localObjVal
println("grbObjVal = $grbObjVal, grbGap = $grbGap, localObjVal = $localObjVal, localGap = $localGap, ObjBnd = $ObjBnd")
