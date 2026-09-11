module Mono
import ..Settings, JuMP, Gurobi, Random

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
_η(rng) = rand(rng, 0.93:1e-4:0.97)
function Model(
    genv, S, T, Δtˈ1h_ratio,
    Ratea, F, CaD,
    EVnode, EVLmax, EVL01, EVEmax, EVEini,
    Bxini, Enode, EEmax, EEini, EPminDiv,
    WPmax,Wscnvec,Wnode,
    LPmax,LmuTup,Ltype,Lnode,
    GPmax,GPref,Gnode
)
    Qyvec, rngv = [JuMP.AffExpr() for s=1:S], map(Random.Xoshiro, 1:S)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variables(m, begin # a superset of the complicating variables
        bES[e=eachindex(EEmax), t=0:T, u=(0,1), s=1:S], Bin
        bx[a=eachindex(EVLmax), t=-1:T, s=1:S], Bin
        0 <= eS[e=eachindex(EEmax), t=-1:T, s=1:S] <= EEmax[e]
        0 <= eV[a=eachindex(EVEmax), t=-1:T, s=1:S] <= (1-max(t,0)/10T)EVEmax[a]
    end)
    
    # Here is the NA constr
    JuMP.@variables(m, begin # common aux variables
        bES_com[e=eachindex(EEmax), t=(0,), u=(0,1), s=(0,)]
        bx_com[a=eachindex(EVLmax), t=(0,), s=(0,)]
        eS_com[e=eachindex(EEmax), t=(0,), s=(0,)]
        eV_com[a=eachindex(EVEmax), t=(0,), s=(0,)]
    end); JuMP.@constraints(m, begin # at t=0, add NA constrs
        [e=eachindex(EEmax), u=(0,1), s=1:S], bES[e,0,u,s] == bES_com[e,0,u,0]
        [a=eachindex(EVLmax), s=1:S], bx[a,0,s] == bx_com[a,0,0]
        [e=eachindex(EEmax), s=1:S], eS[e,0,s] == eS_com[e,0,0]
        [a=eachindex(EVEmax), s=1:S], eV[a,0,s] == eV_com[a,0,0]
    end)

    JuMP.@variables(m, begin # Station variables
        0 <= bu[a=eachindex(EVLmax), t=0:T, s=1:S] <= 1
        0 <= bv[a=eachindex(EVLmax), t=0:T, s=1:S] <= 1
        0 <= pDR[a=eachindex(EVLmax), t=0:T, u=(0,1), s=1:S]
        pEV[a=eachindex(EVLmax), t=0:T, s=1:S]
    end)
    for s=1:S, (a,Pl)=enumerate(EVLmax) # Station inner
        JuMP.set_lower_bound(eV[a,T,s], max(.2 * EVEmax[a], .4 * EVEini[a])) # End-of-horizon level
        JuMP.unset_binary(bx[a,-1,s]); JuMP.fix(bx[a,-1,s], Bxini[a]; force=true) # initial `x` variable
        JuMP.fix(eV[a,-1,s], EVEini[a]; force=true)
        L01vec, ηc = EVL01[a], _η(rngv[s])
        for t=0:T
            JuMP.@constraints(m, begin
                bu[a,t,s]-bv[a,t,s] == bx[a,t,s]-bx[a,t-1,s]
                bx[a,t,s] ≤ 1-bv[a,t,s] # min down time is 1
                sum(bu[a,i,s] for i=range(t;step=-1,length=4+a) if i≥0) ≤ bx[a,t,s] # min up time is UT
                0.1bx[a,t,s]Pl ≤ pEV[a,t,s] # minimum charging power (timeless)
                pEV[a,t,s] ≤ 2.5bx[a,t,s]Pl # maximum charging power (timeless)
                eV[a,t,s]-eV[a,t-1,s] == ηc*(pEV[a,t,s]-L01vec[t+1]Pl)Δtˈ1h_ratio - pDR[a,t,1,s] + pDR[a,t,0,s]
            end)
        end
    end

    JuMP.@variable(m, pES[e=eachindex(EEmax), t=0:T, u=(0,1), s=1:S])
    C01 = 1.5:0.1:6.5 # degradation cost
    for s=1:S, (e,Tp)=enumerate(EPminDiv) # ES unit - inner constrs
        c01 = rand(rngv[s], C01); for t=0:T, u=(0,1)
            JuMP.add_to_expression!(Qyvec[s], c01, pES[e,t,u,s]) # degrade cost
        end
        JuMP.@constraint(m, [t=0:T], bES[e,t,0,s] + bES[e,t,1,s] <= true) # only 3 states
        JuMP.@constraint(m, [t=0:T, u=(0,1)], (EEmax[e]/Tp[u+1])bES[e,t,u,s] <= pES[e,t,u,s]) # Min_power
        JuMP.@constraint(m, [t=0:T, u=(0,1)], pES[e,t,u,s] <= 2.0EEmax[e]bES[e,t,u,s]) # Max_power
        JuMP.set_lower_bound(eS[e,T,s], max(0.2 * EEmax[e], 0.5 * EEini[e]))
        JuMP.fix(eS[e,-1,s], EEini[e]; force=true)
        ηc, ηd = _η(rngv[s]), _η(rngv[s])
        JuMP.@constraint(m, [t=0:T], eS[e,t,s]-eS[e,t-1,s] == (pES[e,t,1,s]ηc-pES[e,t,0,s]/ηd)Δtˈ1h_ratio)
    end

    pfe = Dict((b,t,s) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T for s=1:S)
    gpe = Dict((t,s) => JuMP.AffExpr() for t=0:T for s=1:S)
    for s=1:S, (e,node)=enumerate(Enode), t=0:T
        JuMP.add_to_expression!(gpe[t,s], pES[e,t,0,s])
        JuMP.add_to_expression!(gpe[t,s], -1., pES[e,t,1,s])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t,s],  Fnl, pES[e,t,0,s])
            JuMP.add_to_expression!(pfe[b,t,s], -Fnl, pES[e,t,1,s])
        end
    end
    for s=1:S, (a,node)=enumerate(EVnode), t=0:T
        JuMP.add_to_expression!(gpe[t,s], -1., pEV[a,t,s])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t,s], -Fnl, pEV[a,t,s])
        end
    end
    Cdr = 100:.01:500
    for s=1:S, p=pDR[:, :, :, s] JuMP.add_to_expression!(Qyvec[s], rand(rngv[s], Cdr), p) end

    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T, s=1:S])
    Cϖ = 150:0.01:250
    for s=1:S, p=ϖ[:, :, s] JuMP.add_to_expression!(Qyvec[s], rand(rngv[s], Cϖ), p) end
    for s=1:S, (w,node)=enumerate(Wnode) # wind - to system
        Pmax, Wmat = WPmax[w], Wscnvec[s] # WmuMat
        for t=0:T
            Pw = Pmax * Wmat[1+t, w]
            JuMP.set_upper_bound(ϖ[w,t,s], Pw)
            JuMP.add_to_expression!(gpe[t,s], -1., ϖ[w,t,s])
            JuMP.add_to_expression!(gpe[t,s], Pw)
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t,s], -Fnl, ϖ[w,t,s])
                JuMP.add_to_expression!(pfe[b,t,s],  Fnl * Pw)
            end
        end
    end

    JuMP.@variable(m, 0 <= ζ[l=eachindex(LPmax), t=0:T, s=1:S])
    Cζ = 1000:.1:1800
    for s=1:S, p=ζ[:, :, s] JuMP.add_to_expression!(Qyvec[s], rand(rngv[s], Cζ), p) end
    for s=1:S, (l,node)=enumerate(Lnode) # load - to system
        LdMax, L01Curve = LPmax[l], LmuTup[Ltype[l]]
        for t=0:T
            Pl = LdMax * L01Curve[1+t]
            JuMP.set_upper_bound(ζ[l,t,s], Pl)
            JuMP.add_to_expression!(gpe[t,s], ζ[l,t,s])
            JuMP.add_to_expression!(gpe[t,s], -Pl)
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t,s], Fnl, ζ[l,t,s])
                JuMP.add_to_expression!(pfe[b,t,s], -Fnl * Pl)
            end
        end
    end

    JuMP.@variable(m, 0 <= r[g=eachindex(GPmax), t=0:T, u=(0,1), s=1:S])
    Cr = 7:.1:37
    for s=1:S, p=r[:, :, :, s] JuMP.add_to_expression!(Qyvec[s], rand(rngv[s], Cr), p) end
    for s=1:S, (g,node)=enumerate(Gnode) # Generator Redispatch
        Pmax = GPmax[g]
        for t=0:T
            Pba = Pmax * GPref[t+1,g]
            JuMP.set_upper_bound(r[g,t,0,s], Pba)
            JuMP.set_upper_bound(r[g,t,1,s], Pmax - Pba)
            JuMP.add_to_expression!(gpe[t,s], Pba)
            JuMP.add_to_expression!(gpe[t,s], -1., r[g,t,0,s])
            JuMP.add_to_expression!(gpe[t,s], r[g,t,1,s])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t,s], Fnl * Pba)
                JuMP.add_to_expression!(pfe[b,t,s], -Fnl, r[g,t,0,s])
                JuMP.add_to_expression!(pfe[b,t,s], Fnl, r[g,t,1,s])
            end
        end
    end

    JuMP.@constraint(m, [t=0:T, s=1:S], gpe[t,s] == 0)
    JuMP.@variable(m, 0 <= pϵ[b=eachindex(Ratea), t=0:T, s=1:S])
    for s=1:S, (b,Lim)=enumerate(Ratea)
        d = CaD["branch"][string(b)]
        (d["f_bus"] ∈ Wnode || d["t_bus"] ∈ Wnode) && continue
        JuMP.@constraint(m, [t=0:T],  pfe[b,t,s] <= Lim + pϵ[b,t,s])
        JuMP.@constraint(m, [t=0:T], -Lim -pϵ[b,t,s] <= pfe[b,t,s])
    end
    Cϵ = 9999:.1:11000
    for s=1:S, p=pϵ JuMP.add_to_expression!(Qyvec[s], rand(rngv[s], Cϵ), p) end

    JuMP.@objective(m, Min, sum(Qyvec)/S)
    m
end

end
