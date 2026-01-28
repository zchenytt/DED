module Balance # ✅
import ..WindCur, JuMP

_Wind(t, s, GD) = sum(WindCur.UB(GD,t,s,z,g) for z=GD.Zone for g=eachindex(GD.N[z]))
_Load(t, GD) = sum(GD.P[z][t,g] for z=GD.Zone for g=eachindex(GD.N[z]))
_Const(t, s, Wind, Load, Demand) = _Wind(t, s, Wind) - _Load(t, Load) - _Load(t, Demand)

function _gen(m, t, s, Wind, Demand, p, p0, ϖ, ζ, δ) # for the system at (time=t, scene=s)
    # only `p` and `p0` are from 1st-stage
    e = JuMP.@expression(m, p[t, s] + p0[t, s] +
        sum(ζ[t,s,z,g,0]-ζ[t,s,z,g,1] for z=Demand.Zone for g=eachindex(Demand.N[z])) -
        sum(ϖ[t,s,z,g] for z=Wind.Zone for g=eachindex(Wind.N[z]))
    )
    isnothing(δ) && return e
    e + δ[t, s]
end

add(m, t, s, Wind, Demand, Load, p, p0, ϖ, ζ) = JuMP.@constraint(m, # 1st-stage
    [t=t, s=s], # any instant, any scene
    _gen(m, t, s, Wind, Demand, p, p0, ϖ, ζ, nothing) == -_Const(t, s, Wind, Load, Demand)
)

function add(m, t, s, Wind, Demand, Load, p, p0, ϖ, ζ, δ)
    JuMP.@variable(m, 0 <= a[t=t, s=s, u=(0, 1)] <= 0) # delete_upper_bound(a[t, s, u])
    JuMP.@constraint(m,
        [t=t, s=s], # any instant, any scene
        a[t,s,1] - a[t,s,0] +
        _gen(m, t, s, Wind, Demand, p, p0, ϖ, ζ, δ) == -_Const(t, s, Wind, Load, Demand)
    )
end

end

