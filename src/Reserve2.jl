module Reserve2
import JuMP
# ✅
# There's only 1 type of constr, but split into an initial one and a set of later ones.
# only `r0` is from 1st-stage, defined in Reserve1.jl
_add(m, ramp, UB) = JuMP.@constraint(m, ramp ≤ UB) # ✅
function _add(m, ramp, RampUp, RampDn) # ✅
    _add(m,  ramp, RampUp)
    _add(m, -ramp, RampDn)
end
#=▶=# function add(m, t::Int, s, GD, r0, δ) # `t` is simply the beginning index of planning
    Ramp, t = GD.Ramp, t+1
    for s=s, z=GD.Zone, g=eachindex(GD.N[z]) # each ReserveGenerator amid each scene
        ramp = r0[t,s,z,g] + δ[t,s,z,g]
        _add(m, ramp, Ramp[1][z][g], Ramp[0][z][g])
    end
end
#=▶=# function add(m, t::Int, T::Int, s, GD, r0, δ) # t2T is simply the planning horizon
    Ramp = GD.Ramp
    for t=t+2:T, s=s, z=GD.Zone, g=eachindex(GD.N[z])
        ramp = JuMP.@expression(m, r0[t,s,z,g] + δ[t,s,z,g] - δ[t-1,s,z,g])
        _add(m, ramp, Ramp[1][z][g], Ramp[0][z][g])
    end
end

end

