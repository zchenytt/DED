module Settings
import JuMP, Gurobi

opt_and_ter(m) = (JuMP.optimize!(m); JuMP.termination_status(m))
# The "Crossover=0" option is observed to be too vague (e.g. terminate with a 3% rGap) for certain cases, thus abandoned
# Since Crossover is needed, "Method=2" won't be faster (supposedly)
function reset_param_gap(m)
    JuMP.set_attribute(m, "MIPGap", 1e-4)
    JuMP.set_attribute(m, "MIPGapAbs", 1e-10)
end
const C = Dict{String, Any}("OutputFlag" => 0, "Threads" => 1, "MIPGap" => 0, "MIPGapAbs" => 0)
reset_gurobi_seed(m) = JuMP.set_attribute(m, "Seed", rand(0:2000000000))
function solve_many_times(m, f_where)
    pri = JuMP.INFEASIBLE_POINT
    for k = 0:4
        JuMP.optimize!(m) # you need to set TimeLimit in advance, before calling this function
        pri = JuMP.primal_status(m)
        pri === JuMP.FEASIBLE_POINT && break
        @ccall(printf("%s> unsolved at k=%d\n"::Cstring; f_where::Cstring, k::Cint)::Cint)
        reset_gurobi_seed(m)
    end
    pri === JuMP.FEASIBLE_POINT || error("$f_where: $pri")
end

# The `RC` variable attribute can be used in classic Benders cut generation:
# where the subproblem HAS TO be a _Min_ program.
# Under this restriction, the `RC` attribute of `x` can be seen as `dual(FixRef(x))`

_gcc(o, x) = Gurobi.c_column(o, JuMP.index(x))
function setxdblattrelement(m, x, str, v)
    o = m.moi_backend
    Gurobi.GRBsetdblattrelement(o, str, _gcc(o, x), v)
end
function getxdblattrelement(m, x, str)
    o, r = m.moi_backend, m[:r]
    Gurobi.GRBgetdblattrelement(o, str, _gcc(o, x), r)
    r.x
end
function getmodeldblattr(m, str)
    r = m[:r]
    Gurobi.GRBgetdblattr(m.moi_backend, str, r)
    r.x
end

# Benchmark shows that both `@build_constraint` and `add_constraint` are at the best performance

# API with Expression or Dict should be used judiciously, since oftentimes there
# is better design, e.g. introducing a new variable node (with a == constraint)


# [dual decomposition] used to add `θ - c'β ≤ cn` cuts (signs of θ and cn are regular), `n == length(x/c)`
# [Benders: b = 0/-1 : Feas/MinC] `π'x + (b)θ_s ≤ -Cn`
addle(o, n::Cint, x::Vector{Cint}, c::Vector{Cdouble}, cn::Cdouble) = Gurobi.GRBaddconstr(o, n, x, c, Gurobi.GRB_LESS_EQUAL, cn, "")

Env() = Gurobi.Env(C) # generate a _new_ one as defined by `Gurobi.Env`
Model() = JuMP.direct_model(Gurobi.Optimizer(Env()))
Model(tks, v) = for i=eachindex(tks) setindex!(tks, Threads.@spawn(setindex!(v, Model(), i)), i) end
function Model(tks)
    v = similar(tks, JuMP.Model)
    Model(tks, v)
    foreach(wait, tks)
    v
end

# [Test functions] (monitor the cpu engagement (e.g. by htop) while testing)
function test_multithreaded_mip_build_and_solve(num_threads, N=80)
    tks = Vector{Task}(undef, num_threads)
    println("test> building empty models...")
    v = Model(tks)
    println("test> filling models...")
    build_test_model(tks, v, N)
    println("test> solving models...")
    solve_test_model(tks, v)
end
function solve_test_model(tks, v)
    for (j, m)=enumerate(v) setindex!(tks, Threads.@spawn(JuMP.optimize!(m)), j) end
    foreach(wait, tks)
end
function build_test_model(tks, v, N)
    for (j, m)=enumerate(v) setindex!(tks, Threads.@spawn(build_test_model(m, N)), j) end
    foreach(wait, tks)
end
function build_test_model(m, N)
    # N indicates how hard the MIP is, e.g. N = 80 is hard
    JuMP.@variable(m, x[1:N], Bin)
    JuMP.@variable(m, y[1:N], Bin)
    JuMP.@variable(m, z[1:N, 1:N] >= 0)
    o = m.moi_backend
    for z=z setoc(o, z, rand(-1:.0017:1)) end
    JuMP.set_objective_sense(m, JuMP.MIN_SENSE)
    # add classic BQP cuts
    JuMP.@constraint(m, [i=1:N, j=1:N], z[i, j] <= x[i])
    JuMP.@constraint(m, [i=1:N, j=1:N], z[i, j] <= y[j])
    JuMP.@constraint(m, [i=1:N, j=1:N], z[i, j] >= x[i] + y[j] - 1)
end

printinfo() = (th = map(Threads.nthreads, (:default, :interactive)); println("Settings> Threads=$th"))
end