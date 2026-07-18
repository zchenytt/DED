module Models
import ..Settings, JuMP

function model_and_solve(S, tLeap, T, F, Ratea,
    Wnode, WPmax, Wscnvec,
    LmuTup, Ltype, LPmax, Lnode, LCshd,
    Gnode, GPmax, GClin,
    Enode, EPmax, EEmax, EEini, EClin,
    KCres, KRmp, WCcur,
    Agi, Agp
)
    ge = Settings.Env()
    m = Settings.Model(ge) # monolithic model
    # Linking variable (at t=0, and scene-less)
    JuMP.@variable(m, p0[g=eachindex(Agi)])    # [t=0]
    JuMP.@variable(m, en0[e=eachindex(EEmax)]) # [t=0]
    for (e,Emax)=enumerate(EEmax) # This may be redundant in monolithic solving, but not when decomposing
        Em1, Pma = EEini[e], EPmax[e]
        JuMP.set_lower_bound(en0[e], max(Em1 - Pma, 0.))
        JuMP.set_upper_bound(en0[e], min(Em1 + Pma, Emax))
    end
    for (g,i)=enumerate(Agi)
        Pmax, Pm1 = GPmax[i], Agp[g]
        RmpLim = KRmp * Pmax
        JuMP.set_lower_bound(p0[g], max(Pm1 - RmpLim, 0.))
        JuMP.set_upper_bound(p0[g], min(Pm1 + RmpLim, Pmax))
    end
    JuMP.@variable(m, 0 <= en[e=eachindex(EEmax), t=1:T, s=1:S] <= EEmax[e])
    JuMP.@variable(m, 0 <= pe[e=eachindex(EPmax), t=0:T, u=(0,1), s=1:S] <= EPmax[e])
    JuMP.@variable(m, 0 <=  p[g=eachindex(Agi), t=1:T, s=1:S]) # power output
    JuMP.@variable(m, 0 <=  r[g=eachindex(Agi), t=0:T, u=(0,1), s=1:S]) # Up/Dn reserve 

    for s=1:S, e=eachindex(EPmax) # ES unit - inner dynamics
        local t
        t=0; JuMP.@constraint(m, EEini[e] + 0.95pe[e,t,1,s] - pe[e,t,0,s]/0.95 == en0[e])
        t=1; JuMP.@constraint(m, en0[e]   + 0.95pe[e,t,1,s] - pe[e,t,0,s]/0.95 == en[e,t,s])
        JuMP.@constraint(m, [t=2:T], en[e,t-1,s] + 0.95pe[e,t,1,s] - pe[e,t,0,s]/0.95 == en[e,t,s])
        JuMP.set_lower_bound(en[e,T,s], EEini[e])
    end
    for s=1:S, (g,i)=enumerate(Agi) # generator inner dynamics
        Pmax = GPmax[i]; PRmp = KRmp * Pmax
        JuMP.@constraint(m, p0[g] - r[g,0,0,s] >= 0)
        JuMP.@constraint(m, p0[g] + r[g,0,1,s] <= Pmax)
        JuMP.@constraint(m, [t=1:T], p[g,t,s] - r[g,t,0,s] >= 0)
        JuMP.@constraint(m, [t=1:T], p[g,t,s] + r[g,t,1,s] <= Pmax)
        JuMP.@constraint(m, p[g,1,s] <= p0[g] + PRmp)
        JuMP.@constraint(m, p[g,1,s] >= p0[g] - PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t,s] <= p[g,t-1,s] + PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t,s] >= p[g,t-1,s] - PRmp)
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,1,s]) >= 4) # system-wide Up reserve demand
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,0,s]) >= 2) # system-wide Dn reserve demand
    end

    Obj2 = [JuMP.AffExpr() for s=1:S]; # System-wide constructs
    JuMP.@variable(m, -Ratea[b] <= pf[b=eachindex(Ratea), t=0:T, s=1:S] <= Ratea[b]);
    pfe = Dict((b,t,s) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T for s=1:S);
    gpe = Dict((t,s)   => JuMP.AffExpr()                        for t=0:T for s=1:S);
    for s=1:S, (e,node)=enumerate(Enode), t=0:T # ES unit - to system
        JuMP.add_to_expression!(Obj2[s], EClin, pe[e,t,0,s])
        JuMP.add_to_expression!(Obj2[s], EClin, pe[e,t,1,s])
        JuMP.add_to_expression!(gpe[t,s], pe[e,t,0,s])
        JuMP.add_to_expression!(gpe[t,s], -1., pe[e,t,1,s])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t,s],  Fnl, pe[e,t,0,s])
            JuMP.add_to_expression!(pfe[b,t,s], -Fnl, pe[e,t,1,s])
        end
    end
    for s=1:S, (g,i)=enumerate(Agi) # generator - to system
        gC, node = GClin[i], Gnode[i]; rC = KCres * gC
        for t=0:T JuMP.add_to_expression!(Obj2[s], rC, r[g,t,0,s]) end
        for t=0:T JuMP.add_to_expression!(Obj2[s], rC, r[g,t,1,s]) end
        local t
        t=0; JuMP.add_to_expression!(Obj2[s], gC, p0[g])
        for t=1:T JuMP.add_to_expression!(Obj2[s], gC, p[g,t,s]) end
        t=0; JuMP.add_to_expression!(gpe[t,s], p0[g])
        for t=1:T JuMP.add_to_expression!(gpe[t,s], p[g,t,s]) end
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            t=0; JuMP.add_to_expression!(pfe[b,t,s], Fnl, p0[g])
            for t=1:T JuMP.add_to_expression!(pfe[b,t,s], Fnl, p[g,t,s]) end
        end
    end
    JuMP.@variable(m, 0 <= ζ[l=eachindex(Lnode), t=0:T, s=1:S]) # load shedding
    for (l,node)=enumerate(Lnode) # load - to system
        ty = Ltype[l]
        for t=0:T
            Pmax, l01 = LPmax[l], LmuTup[ty][1+(tLeap)t]
            Pl = Pmax * l01
            for s=1:S
                JuMP.set_upper_bound(ζ[l,t,s], Pl)
                JuMP.add_to_expression!(gpe[t,s], -Pl)
                JuMP.add_to_expression!(gpe[t,s], ζ[l,t,s])
                JuMP.add_to_expression!(Obj2[s], LCshd, ζ[l,t,s])
                for b=eachindex(Ratea)
                    Fnl = F[b, node]
                    JuMP.add_to_expression!(pfe[b,t,s], -Fnl * Pl)
                    JuMP.add_to_expression!(pfe[b,t,s], Fnl, ζ[l,t,s])
                end
            end
        end
    end
    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T, s=1:S]);
    for (s,Wmat)=enumerate(Wscnvec), (w,node)=enumerate(Wnode)
        Pmax = WPmax[w]
        for t=0:T
            w01 = Wmat[1+(tLeap)t, w]
            Pw = Pmax * w01
            JuMP.set_upper_bound(ϖ[w,t,s], Pw)
            JuMP.add_to_expression!(gpe[t,s], Pw)
            JuMP.add_to_expression!(gpe[t,s], -1., ϖ[w,t,s])
            JuMP.add_to_expression!(Obj2[s], WCcur, ϖ[w,t,s])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t,s],  Fnl * Pw)
                JuMP.add_to_expression!(pfe[b,t,s], -Fnl, ϖ[w,t,s])
            end
        end
    end

    JuMP.@constraint(m, [t=0:T, s=1:S], gpe[t,s] == 0)
    JuMP.@constraint(m, [b=eachindex(Ratea), t=0:T, s=1:S], pfe[b,t,s] == pf[b,t,s])
    JuMP.@objective(m, Min, sum(Obj2)/S)
    
    # JuMP.set_attribute(m, "PreDual", 0)
    # JuMP.set_attribute(m, "Presolve", 0)
    # JuMP.set_attribute(m, "Method", 0)
    JuMP.set_attribute(m, "Threads", 2) # JuMP.set_time_limit_sec(m, 3600.0)
    JuMP.unset_silent(m)

    m, Obj2, p0, en0, en, pe, p, r, ϖ, ζ
end

end
