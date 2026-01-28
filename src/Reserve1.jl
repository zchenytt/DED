module Reserve1
import JuMP
# ✅
function _add(m, p0, δ0up, δ0dn, UB)
    JuMP.@constraint(m, p0 + δ0up ≤ UB)
    JuMP.@constraint(m, 0.0 ≤ p0 - δ0dn) # TODO: add a nonzero PMin
end
#=▶=# function add(m, t, s, GD, p0, δ0)
    PMax = GD.PMax
    for t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z])
        _add(m, p0[t,s,z,g], δ0[t,s,z,g,1], δ0[t,s,z,g,0], PMax[z][g]) 
    end
end

# This 1st-stage ramping constr motivates introducing `r0`
function _add(m, r0t, δ0tUp, δ0tDn, RampUp, RampDn)
    JuMP.@constraint(m, r0t + δ0tUp ≤ RampUp)
    JuMP.@constraint(m, δ0tDn - r0t ≤ RampDn)
end
#=▶=# function add_r0(m, t::Int, T::Int, s, GD, p0) # totally decided by `p0`, so ∈ 1st-stage, too
    # indexed by the latter time index, aligning with the <1st-stage> ramp constrs
    r0 = JuMP.@variable(m, [t=t+1:T, s=s, z=GD.Zone, g=eachindex(GD.N[z])])
    JuMP.@constraint(m,    [t=t+1:T, s=s, z=GD.Zone, g=eachindex(GD.N[z])],
        r0[t,s,z,g] == p0[t,s,z,g] - p0[t-1,s,z,g] # Definition
    )
    r0
end
#=▶=# function add(m, t::Int, T::Int, s, GD, r0, δ0) # this t2T is simply the planning horizon
    Ramp = GD.Ramp
    for t=t+1:T, s=s, z=GD.Zone, g=eachindex(GD.N[z])
        _add(m, r0[t,s,z,g], δ0[t,s,z,g,1], δ0[t,s,z,g,0], Ramp[1][z][g], Ramp[0][z][g])
    end
end

end
