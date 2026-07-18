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

# Below are test results
julia> JuMP.optimize!(m)
Gurobi Optimizer version 13.0.1 build v13.0.1rc0 (linux64gpu - "Ubuntu 24.04.4 LTS")

CPU model: AMD EPYC 7763 64-Core Processor, instruction set [SSE2|AVX|AVX2]
Thread count: 128 physical cores, 256 logical processors, using up to 2 threads

GPU model: NVIDIA RTX A6000, CUDA compute version 8.6, NVIDIA driver compatible with CUDA version 13

Non-default parameters:
Threads  2

Optimize a model with 15485000 rows, 21000012 columns and 1277800000 nonzeros (Min)
Model fingerprint: 0xfc8edb72
Model has 8790010 linear objective coefficients
Coefficient statistics:
  Matrix range     [2e-06, 1e+00]
  Objective range  [3e-05, 3e+00]
  Bounds range     [1e-02, 7e+01]
  RHS range        [2e-07, 5e+01]

Presolve removed 4884604 rows and 6191969 columns (presolve time = 1113s)...
Presolve time: 1304.50s
Presolved: 10600396 rows, 14808043 columns, 899710611 nonzeros

Concurrent LP optimizer: dual simplex and barrier
Showing barrier log only...

Elapsed ordering time = 14s
Ordering time: 102.31s

Barrier statistics:
 Dense cols : 12
 AA' NZ     : 9.704e+08
 Factor NZ  : 1.137e+09 (roughly 20.0 GB of memory)
 Factor Ops : 1.598e+11 (roughly 30 seconds per iteration)
 Threads    : 1

                  Objective                Residual
Iter       Primal          Dual         Primal    Dual     Compl     Time
   0   2.41749670e+07 -2.79482073e+08  1.62e+04 6.21e-01  9.68e+03  1788s
  79   1.94861956e+04  1.94861954e+04  5.17e-07 6.19e-07  3.25e-10  6659s

Barrier solved model in 79 iterations and 6658.93 seconds (15183.61 work units)
Optimal objective 1.94861956e+04

Crossover log...
 9107016 variables added to crossover basis                     6754s
 1291956 DPushes remaining with DInf 0.0000000e+00              6766s
       0 DPushes remaining with DInf 0.0000000e+00              7237s
 3902213 PPushes remaining with PInf 7.5016081e-02              7242s
       0 PPushes remaining with PInf 0.0000000e+00              8025s
  Push phase complete: Pinf 0.0000000e+00, Dinf 4.3775332e+01   8027s

Iteration    Objective       Primal Inf.    Dual Inf.      Time
 5192890    1.9486196e+04   0.000000e+00   4.377533e+01   8056s
Crossover time: 1400.09 seconds (884.89 work units)

Solved with barrier
 5192929    1.9486195e+04   0.000000e+00   0.000000e+00   8141s

Solved in 5192929 iterations and 8140.91 seconds (16106.80 work units)
Optimal objective  1.948619542e+04

User-callback calls 457887, time in user-callback 0.39 sec

julia> for s=1:S
           all(x -> abs(x) < 1e-6, JuMP.value.(pe[:, :, 1, s]) .* JuMP.value.(pe[:, :, 0, s])) || @error("s = $s")
       end
┌ Error: s = 833
└ @ Main REPL[37]:2
...
┌ Error: s = 4948
└ @ Main REPL[37]:2

julia> JuMP.value.(pe[:, :, 0, 4948])
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 0:12
And data, a 2×13 Matrix{Float64}:
 0.0  0.0        0.0       0.0  0.0  0.0  2.69781  2.93631  1.76736  2.44776  5.31749  0.540776  0.654784
 0.0  0.0550276  0.603883  0.0  0.0  0.0  0.0      0.0      0.0      0.0      0.0      0.0       0.95

julia> JuMP.value.(pe[:, :, 1, 4948])
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 0:12
And data, a 2×13 Matrix{Float64}:
 3.90112  6.30838  4.76794   2.96463  0.187888  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0      1.1136   0.669123  0.0      0.0       0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> extrema(JuMP.value.(Obj2))
(16512.882718150362, 23213.51546543862)

julia> import Statistics; Statistics.mean(JuMP.value.(Obj2))
19486.195419125313

julia> JuMP.value.(p0)'
1×10 adjoint(::Vector{Float64}) with eltype Float64:
 1.89375  4.85  0.125  2.23  1.54  1.95  2.89325  3.185  6.53  0.229

julia> JuMP.value.(en0)'
1×2 adjoint(::Vector{Float64}) with eltype Float64:
 4.70606  34.0

