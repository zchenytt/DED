module Get
import JuMP, Gurobi

# Retrieve the dual solution of copy constrs
_π(#=dual ∘ FixRef=# dual, #=inn[s]=# m) = (
    p  = dual.(m[:p ]),
    p0 = dual.(m[:p0]),
    r0 = dual.(m[:r0]),
    δ0 = dual.(m[:δ0])
)
fix(m, X, properties) = for p=properties foreach(JuMP.fix, m[p], X[p]) end
function innerstats!(exp, NegObj, m, x_link, X, properties)
    for p=properties, (zCopy, x, X) = zip(m[p], x_link[p], X[p])
        pi = JuMP.dual(JuMP.FixRef(zCopy))
        NegObj += (pi)X
        JuMP.add_to_expression!(exp, pi, x)
    end
    NegObj # `exp` is returned implicitly
end

######################################
# Retrieve trial solution && Exporting
######################################
linkingSysPower(p, t) = view(JuMP.value.(p), t, [1])
function trial_solution_nt(mst, θ, #=t+1:T=# t, _p, _p0, _r0, _δ0)
    p = linkingSysPower(_p, t)
    p0 = linkingSysPower(_p0, t)
    r0 = JuMP.value.(_r0)
    δ0 = JuMP.value.(_δ0)
    θ = JuMP.value.(θ)
    lb = JuMP.objective_bound(mst)
    x = (; p, p0, r0, δ0)
    (; lb, θ, x)
end

end
