module Settings
import JuMP, Gurobi

const C = Dict{String, Any}("OutputFlag" => 0, "Threads" => 1, "Method" => 2, "Crossover" => 0)
Env() = Gurobi.Env(C) # generate a _new_ one as defined by `Gurobi.Env`
Model() = JuMP.direct_model(Gurobi.Optimizer(Env()))

function _fill!(v, i)
    v[i] = Model()
    return
end
_undef(N) = Vector{JuMP.Model}(undef, N)
function Model(N) # still faster than 1-thread serial construction, i.e. [Model() for i=1:N]
    v = _undef(N)
    foreach(wait, [Threads.@spawn(_fill!(v, i)) for i=1:N])
    v
end

end
