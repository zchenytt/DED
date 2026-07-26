"""
This pre-model doesn't consider >1 scenarios.

To decide the subset of generators that is at ON state during the planning
and also decide initial conditions (e.g. Energy/Power_(t = -1))
"""
module PreModel

ge = Settings.Env();
m = Settings.Model(ge);
Obj = JuMP.AffExpr();
pfe = Dict((b, t) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T);
gpe = Dict(t => JuMP.AffExpr()                             for t=0:T);
JuMP.@variable(m, 0 <= p[g=eachindex(GPmax), t=0:T] <= GPmax[g]); # It can happen that sum(Wind in sys) < sum(Load in sys), so we need power from generators
JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T]); # It can happen that sum(Wind in sys) alone > sum(Load in sys), so we need wind curtailments
JuMP.@variable(m, 0 <= pe[e=eachindex(EPmax), t=0:T, u=(0,1)] <= EPmax[e]); # discharge_0/charge_1 power
JuMP.@variable(m, 0 <= en[e=eachindex(EEmax), t = -1:T] <= EEmax[e]);
EEini = [en[e,-1] for e=eachindex(EEmax)]; # Follows ES unit - inner dynamics
JuMP.@constraint(m, [e=eachindex(EEmax), t=0:T], en[e,t-1] + 0.95pe[e,t,1] - pe[e,t,0]/0.95 == en[e,t])
JuMP.@constraint(m, [e=eachindex(EEmax)], en[e,T] >= en[e,-1])
for (e,node)=enumerate(Enode), t=0:T # ES unit - outer behavior
    JuMP.add_to_expression!(Obj, EClin, pe[e,t,0])
    JuMP.add_to_expression!(Obj, EClin, pe[e,t,1])
    JuMP.add_to_expression!(gpe[t], pe[e,t,0])
    JuMP.add_to_expression!(gpe[t], -1., pe[e,t,1])
    for b=eachindex(Ratea)
        Fnl = F[b, node]
        JuMP.add_to_expression!(pfe[b,t],  Fnl, pe[e,t,0])
        JuMP.add_to_expression!(pfe[b,t], -Fnl, pe[e,t,1])
    end
end
for (g,node)=enumerate(Gnode)
    Rmp_g = KRmp * GPmax[g]
    JuMP.@constraint(m, [t=1:T], p[g,t] - p[g,t-1] <= Rmp_g)
    JuMP.@constraint(m, [t=1:T], p[g,t-1] - p[g,t] <= Rmp_g)
    for t=0:T
        JuMP.add_to_expression!(Obj, GClin[g], p[g,t])
        JuMP.add_to_expression!(gpe[t],    p[g,t])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t], Fnl, p[g,t]) # we haven't introduced p_injection yet, for simplicity
        end
    end
end
for (w,node)=enumerate(Wnode), t=0:T
    Pmax, w01 = WPmax[w], WmuMat[1+(tLeap)t, w]
    Pw = Pmax * w01
    JuMP.set_upper_bound(ϖ[w,t], Pw)
    JuMP.add_to_expression!(gpe[t], Pw)
    JuMP.add_to_expression!(gpe[t], -1., ϖ[w,t])
    JuMP.add_to_expression!(Obj, WCcur, ϖ[w,t])
    for b=eachindex(Ratea)
        Fnl = F[b, node]
        JuMP.add_to_expression!(pfe[b, t], Fnl * Pw)
        JuMP.add_to_expression!(pfe[b, t], -Fnl, ϖ[w,t])
    end
end
for (l,node)=enumerate(Lnode), t=0:T
    Pmax, l01 = LPmax[l], LmuTup[Ltype[l]][1+(tLeap)t]
    Pl = Pmax * l01
    JuMP.add_to_expression!(gpe[t], -Pl)
    for b=eachindex(Ratea)
        Fnl = F[b, node]
        JuMP.add_to_expression!(pfe[b, t], -Fnl * Pl)
    end
end
JuMP.@variable(m, -Ratea[b] <= pf[b=eachindex(Ratea), t=0:T] <= Ratea[b]);
JuMP.@constraint(m, [b=eachindex(Ratea), t=0:T], pfe[b, t] == pf[b, t]);
JuMP.@constraint(m, [t=0:T], gpe[t] == 0);
JuMP.@objective(m, Min, Obj);
JuMP.unset_silent(m)
JuMP.optimize!(m)

end
