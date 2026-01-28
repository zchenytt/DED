module Costs
import JuMP

function add_θ_set_obj(m, S, cˈx)
    θ = JuMP.@variable(m, [1:S])
    JuMP.@objective(m, Min, cˈx + 1/S * sum(θ))
    θ
end

_c(D, z, g, u) = (1-u + (u)D.Mul[z][g])D.Slope[z][g]
_p(m, t, s, GD, p) = JuMP.@expression(m,sum(
    _c(GD.Cost.p, z,g,u)p[t,s,z,g,u]
    for t=t for s=s for z=GD.Zone for g=eachindex(GD.N[z]) for u=(0,1)
))
_δ(m, t, s, GD, p) = JuMP.@expression(m,sum(
    _c(GD.Cost.p, z,g,1)p[t,s,z,g,1] # only UpWard Redispatch are priced
    for t=t for s=s for z=GD.Zone for g=eachindex(GD.N[z])
))

_ζ(m, t, s, GD, p) = JuMP.@expression(m,sum(
    GD.Cost.p[u][z][t, g]p[t,s,z,g,u] # DR incentives are paid
    for t=t for s=s for z=GD.Zone for g=eachindex(GD.N[z]) for u=(0,1)
))
_δ0(m, t, s, GD, p) = JuMP.@expression(m,sum(
    GD.Cost.r[u][z][g]p[t,s,z,g,u]
    for t=t for s=s for z=GD.Zone for g=eachindex(GD.N[z]) for u=(0,1)
))

# Used both directly or innerly
_1st_stage(m, t::Int, T::Int, s, Other, Reserve, Demand, p, p0, δ0, ζ) = JuMP.@expression(m,
    _p(m, t:T, s, Other, p) + _p(m, t:T, s, Reserve, p0) +
    _δ0(m, t+1:T, s, Reserve, δ0) +
    _ζ(m, t, s, Demand, ζ)
)
# Multi-period, without any other coefficient
_2nd_stage(m, t::UnitRange{Int}, s, Reserve, Demand, δ, ζ) = _δ(m, t, s, Reserve, δ) + _ζ(m, t, s, Demand, ζ)
# An application instance
add_C1_obj(m, t::Int, T::Int, s, Other, Reserve, Demand, p, p0, δ0, ζ, δ) = JuMP.@objective(
    m, Min,
    _2nd_stage(m, t+1:T, s, Reserve, Demand, δ, ζ) +
    _1st_stage(m, t, T, s, Other, Reserve, Demand, p, p0, δ0, ζ)
)


end

