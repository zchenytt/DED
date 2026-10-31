"Traditional models for online replanning purpose"
module Mono2
import ..Settings, JuMP, Gurobi, Random, ..WindGen

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
_η(rng) = rand(rng, 0.93:1e-4:0.97)
function buildModel(modelrng, # only build once, and then modify-and-resolve during online replanning
    genv, S, T, Δtˈ1h_ratio,
    Ratea, F, CaD,
    EVnode, EVLmax, EVL01, EVEmax, EVEini,
    Bxini, Enode, EEmax, EEini, EPminDiv,
    WPmax,Wscnvec,Wnode,
    LPmax,LmuTup,Ltype,Lnode,
    GPmax,GPref,Gnode)

    QY = Dict(t => NaN for t=0:T) # 🟧used to collect solutions of online replanning
    Qy = Dict((t, s) => JuMP.AffExpr() for s=1:S for t=1:T) # 🟧update
    Qy[0, 1] = JuMP.AffExpr()
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variables(m, begin # a superset of the complicating variables
        bES[e=eachindex(EEmax), t=0:T, u=(0,1), s=1:S], Bin
        bx[a=eachindex(EVLmax), t=-1:T, s=1:S], Bin
        0 <= eS[e=eachindex(EEmax), t=-1:T, s=1:S] <= EEmax[e]
        0 <= eV[a=eachindex(EVEmax), t=-1:T, s=1:S] <= (1-max(t,0)/10T)EVEmax[a]
    end)
    bESX = Dict((e,t,u) => NaN for e=eachindex(EEmax) for t=0:T for u=(0,1))
    bxX = Dict((a,t) => NaN for a=eachindex(EVLmax) for t=0:T)
    buX = copy(bxX)
    eSX = Dict((e,t) => NaN for e=eachindex(EEmax) for t=0:T)
    eVX = Dict((a,t) => NaN for a=eachindex(EVEmax) for t=0:T)
    JuMP.@constraints(m, begin # 🟧update
        [e=eachindex(EEmax), t=0:T, u=(0,1), s=2:S], bES[e,t,u,s] == bES[e,t,u,1]
        [a=eachindex(EVLmax), t=0:T, s=2:S], bx[a,t,s] == bx[a,t,1]
        [e=eachindex(EEmax), t=0:T, s=2:S], eS[e,t,s] == eS[e,t,1]
        [a=eachindex(EVEmax), t=0:T, s=2:S], eV[a,t,s] == eV[a,t,1]
    end) # Here is the NA constr at t=0
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
        L01vec, ηc = EVL01[a], _η(Random.Xoshiro(modelrng + a))
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
        c01 = rand(Random.Xoshiro(modelrng + e), C01)
        for u=(0,1)
            for t=1:T
                JuMP.add_to_expression!(Qy[t, s], c01, pES[e,t,u,s]) # degrade cost
            end
            s==1 && JuMP.add_to_expression!(Qy[0, s], c01, pES[e,0,u,s]) # degrade cost
        end
        JuMP.@constraint(m, [t=0:T], bES[e,t,0,s] + bES[e,t,1,s] <= true) # only 3 states
        JuMP.@constraint(m, [t=0:T, u=(0,1)], (EEmax[e]/Tp[u+1])bES[e,t,u,s] <= pES[e,t,u,s]) # Min_power
        JuMP.@constraint(m, [t=0:T, u=(0,1)], pES[e,t,u,s] <= 2.0EEmax[e]bES[e,t,u,s]) # Max_power
        JuMP.set_lower_bound(eS[e,T,s], max(0.2 * EEmax[e], 0.5 * EEini[e]))
        JuMP.fix(eS[e,-1,s], EEini[e]; force=true)
        ηc = ηd = _η(Random.Xoshiro(modelrng + 2e))
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
    for a=eachindex(EVLmax), u=(0,1), s=1:S
        for t=1:T
            JuMP.add_to_expression!(Qy[t, s], rand(Random.Xoshiro(t+a+u+modelrng), Cdr), pDR[a, t, u, s])
        end
        if s == 1
            JuMP.add_to_expression!(Qy[0, s], rand(Random.Xoshiro(0+a+u+modelrng), Cdr), pDR[a, 0, u, s])
        end
    end

    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T, s=1:S])
    JuMP.@variable(m, pw[w=eachindex(WPmax), t=0:T, s=1:S]) # 🟧 always be fixed to the wind scenes
    JuMP.@constraint(m, ϖ .<= pw)
    Cϖ = 150:0.01:250
    for s=1:S
        for t=1:T
            for p=ϖ[:, t, s] 
                JuMP.add_to_expression!(Qy[t, s], rand(Random.Xoshiro(modelrng+t), Cϖ), p) 
            end
        end
        if s == 1
            for p=ϖ[:, 0, s] 
                JuMP.add_to_expression!(Qy[0, s], rand(Random.Xoshiro(modelrng+0), Cϖ), p) 
            end
        end
    end
    for s=1:S, (w,node)=enumerate(Wnode) # wind - to system
        Pmax, Wmat = WPmax[w], Wscnvec[s] # 🟧🟧
        for t=0:T
            Pw = Pmax * Wmat[1+t, w]
            JuMP.fix(pw[w,t,s], Pw; force=true) # 🟧🟧 to be modified
            JuMP.add_to_expression!(gpe[t,s], -1., ϖ[w,t,s])
            JuMP.add_to_expression!(gpe[t,s], pw[w,t,s])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t,s], -Fnl, ϖ[w,t,s])
                JuMP.add_to_expression!(pfe[b,t,s],  Fnl, pw[w,t,s])
            end
        end
    end

    JuMP.@variable(m, 0 <= ζ[l=eachindex(LPmax), t=0:T, s=1:S])
    Cζ = 1000:.1:1800
    for s=1:S
        for t = 1:T
            for p=ζ[:, t, s]
                JuMP.add_to_expression!(Qy[t,s], rand(Random.Xoshiro(modelrng+t), Cζ), p) 
            end
        end
        if s == 1
            for p=ζ[:, 0, s]
                JuMP.add_to_expression!(Qy[0,s], rand(Random.Xoshiro(modelrng+0), Cζ), p) 
            end
        end
    end
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
    Cr = 7:.1:37 # figure 22 ≈ Some USD/MWh, see Case118.jl file
    for g=eachindex(GPmax), u=(0,1), s=1:S
        for t = 1:T
            JuMP.add_to_expression!(Qy[t,s], rand(Random.Xoshiro(t+g+u+modelrng), Cr), r[g, t, u, s])
        end
        if s == 1
            JuMP.add_to_expression!(Qy[0,s], rand(Random.Xoshiro(0+g+u+modelrng), Cr), r[g, 0, u, s])
        end
    end

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
    for s=1:S
        for t=1:T
            for p=pϵ[:, t, s]
                JuMP.add_to_expression!(Qy[t, s], rand(Random.Xoshiro(modelrng+t), Cϵ), p)
            end
        end
        if s == 1
            for p=pϵ[:, 0, s]
                JuMP.add_to_expression!(Qy[0, s], rand(Random.Xoshiro(modelrng+0), Cϵ), p)
            end
        end
    end
    JuMP.@objective(m, Min, sum(Qy[t,s]/S for s=1:S for t=1:T) + Qy[0,1]) # 🟧
    (; m, Qy, QY, bES, bx, eS, eV, bu, pw, bESX, bxX, eSX, eVX, buX)
end

function online_replanning!(mnt, zmat, Wscnvec, S, T, TimeLimit_per_SP, EEmax, EVLmax, EVEmax, Wnode, WPmax)
    JuMP.unset_silent(mnt.m)
    JuMP.set_time_limit_sec(mnt.m, 15 * 60.)
    t = 0
    JuMP.optimize!(mnt.m) # initial optimization to decide all the Bin_and_state vars
    mnt.QY[t] = JuMP.value(mnt.Qy[0,1]) # record the immediate cost
    for t = 0:T
        for e=eachindex(EEmax), u=(0,1)
            mnt.bESX[e,t,u] = JuMP.value(mnt.bES[e,t,u,1]) 
        end
        for a=eachindex(EVLmax)
            mnt.bxX[a,t] = JuMP.value(mnt.bx[a,t,1])
            mnt.buX[a,t] = JuMP.value(mnt.bu[a,t,1])
        end
        for e=eachindex(EEmax)
            mnt.eSX[e,t] = JuMP.value(mnt.eS[e,t,1])
        end
        for a=eachindex(EVEmax)
            mnt.eVX[a,t] = JuMP.value(mnt.eV[a,t,1])
        end
    end
    for t = 1:T+1
        for e=eachindex(EEmax), u=(0,1)
            JuMP.fix(mnt.bES[e,t-1,u,1], mnt.bESX[e,t-1,u]; force=true)
        end
        for a=eachindex(EVLmax)
            JuMP.fix(mnt.bx[a,t-1,1], mnt.bxX[a,t-1]; force=true)
            JuMP.fix(mnt.bu[a,t-1,1], mnt.buX[a,t-1]; force=true)
        end
        for e=eachindex(EEmax)
            JuMP.fix(mnt.eS[e,t-1,1], mnt.eSX[e,t-1]; force=true)
        end
        for a=eachindex(EVEmax)
            JuMP.fix(mnt.eV[a,t-1,1], mnt.eVX[a,t-1]; force=true)
        end
    end 
    while true
        t += 1 # 🟦 Fold the horizon
        if t > T
            println("Reached the end of the horizon at t=T=$t")
            break
        end
        WindGen.reveal!(zmat, t, Random.Xoshiro(0))
        WindGen.fill_Wscnvec_once!(Wscnvec, zmat, t+1, S)
        for s=1:S, (w,node)=enumerate(Wnode) # wind - to system
            Pmax, Wmat = WPmax[w], Wscnvec[s] # 🟧🟧
            for t=0:T
                Pw = Pmax * Wmat[1+t, w]
                JuMP.fix(mnt.pw[w,t,s], Pw; force=true) # 🟧🟧 to be modified
            end
        end
        mnt.Qy[0,1] = JuMP.@expression(mnt.m, sum(mnt.Qy[t,s]/S for s=1:S)) # store the immediate cost
        JuMP.@objective(mnt.m, Min, sum(mnt.Qy[t,s]/S for s=1:S for t=t+1:T) + mnt.Qy[0,1])
        JuMP.set_time_limit_sec(mnt.m, TimeLimit_per_SP)
        JuMP.optimize!(mnt.m) # initial optimization to decide all the Bin_and_state vars
        mnt.QY[t] = JuMP.value(mnt.Qy[0,1]) # record the immediate cost
    end
end

end
