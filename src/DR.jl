module DR
import JuMP

# u = 0 => cut down load => ub = P/1
# u = 1 => increase load => ub = P/2
UB(P, t, z, g, u) = P[z][t, g]/(u + 1)
function add(model, t, s, GD)
    P = GD.P
    JuMP.@variable(model,
        [t=t, s=s(t), z=GD.Zone, g=eachindex(GD.N[z]), u=(0,1)],
        lower_bound = 0,
        upper_bound = UB(P, t, z, g, u)
    )
end

end

# ζ = DR.add(model, 1:T, _ -> [s], Demand)
