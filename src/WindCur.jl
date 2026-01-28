module WindCur
import JuMP

UB(GD, t, s, z, g) = GD.S[s][t,z]GD.PMax[z][g]

# ϖ = WindCur.add(model, t:T, _ -> [s], Wind)
add(model, t, s, GD) = JuMP.@variable(model,
    [t=t, s=s(t), z=GD.Zone, g=eachindex(GD.N[z])],
    lower_bound = 0,
    upper_bound = UB(GD, t, s, z, g)
)

end

