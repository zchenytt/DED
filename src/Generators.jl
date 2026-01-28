module Generators # ✅
import JuMP

_ub(GD, z, g, u) = (u + (1 - 2*u)GD.Cost.p.Knee[z][g])GD.PMax[z][g]
_addc_has_ub(m, t, s, GD) = JuMP.@variable(m,
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z]), u=(0,1)],
    lower_bound = 0, upper_bound = _ub(GD, z, g, u)
) # [1st-Stage] Generation
_addc_no_ub(m, t, s, GD) = JuMP.@variable(m, 
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z]), u=(0,1)],
    lower_bound = 0 # ub is in `add_δ_le_δ0_constr`
) # [2nd-stage] used in redispatch costs

_cis(_addc, m, t, s, GD) = (
    c = _addc(m, t, s, GD),
    i = _addi(m, t, s, GD),
    s = _adds(m, t, s)
) # Build a NamedTuple
_addi(m, t, s, GD) = JuMP.@variable(m, [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z])])
_adds(m, t, s)     = JuMP.@variable(m, [t=t, s=s]) # power system-level

_c2i_addition(m, t, s, GD, p) = JuMP.@constraint(m,
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z])],
    p.i[t,s,z,g] == p.c[t,s,z,g,1] + p.c[t,s,z,g,0] # TODO: add a nonzero PMin
) # intra-1st-stage
_c2i_subtract(m, t, s, GD, p) = JuMP.@constraint(m,
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z])],
    p.i[t,s,z,g] == p.c[t,s,z,g,1] - p.c[t,s,z,g,0]
) # intra-2nd-stage

# Calling this from another module:
# p/p0 = Generators.add(_addc_has_ub, _c2i_addition, m, t, s, Other/Reserve)
# δ = Generators.add(_addc_no_ub, _c2i_subtract, m, t, s, Reserve)
#=▶=# function add(_addc, _c2i_constr, m, t, s, GD)
    p = _cis(_addc, m, t, s, GD)
    _c2i_constr(m, t, s, GD, p)
    _i2s(m, t, s, GD, p)
    p
end
_i2s(m, t, s, GD, p) = JuMP.@constraint(m, [t=t, s=s],
    p.s[t,s] == sum(p.i[t,s,z,g] for z=GD.Zone for g=eachindex(GD.N[z]))
) # intra 1st/2nd (respectively) stage

# [1st-Stage] Up/Dn Reserve
add_reserve(m, t, s, GD) = JuMP.@variable(m,
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z]), u=(0,1)],
    lower_bound = 0, upper_bound = GD.Ramp[u][z][g]
)

# [2nd-linking-to-1st Stage Constr]
add_δ_le_δ0_constr(m, t, s, GD, δ, δ0) = JuMP.@constraint(m,
    [t=t, s=s, z=GD.Zone, g=eachindex(GD.N[z]), u=(0,1)],
    δ[t,s,z,g,u] ≤ δ0[t,s,z,g,u]
)

end
