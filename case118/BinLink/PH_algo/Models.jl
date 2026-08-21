module Models
import DataStructures.CircularBuffer
import ..Settings, JuMP, Gurobi, Random
const SubMIPTy=@NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, xs::Vector{JuMP.VariableRef}, Qy::JuMP.AffExpr, N::Int64, Xl::Vector{Float64}, Pi::Vector{Float64}, Xl2::Vector{Float64}, Cd::Vector{Float64}}

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
powermat_to_sys(vec, mat, nodevec, T, gpe, pfe, Ratea, F) = for (g,node)=enumerate(nodevec)
    Pmax = vec[g]
    for t=0:T
        Pw = Pmax * mat[1+t, g]
        JuMP.add_to_expression!(gpe[t], Pw)
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t], Fnl * Pw)
        end
    end
end
ES_2_sys(Enode, Ratea, F, T, gpe, pfe, pES) = for (e,node)=enumerate(Enode), t=0:T
    JuMP.add_to_expression!(gpe[t], pES[e,t,0])
    JuMP.add_to_expression!(gpe[t], -1., pES[e,t,1])
    for b=eachindex(Ratea)
        Fnl = F[b, node]
        JuMP.add_to_expression!(pfe[b,t],  Fnl, pES[e,t,0])
        JuMP.add_to_expression!(pfe[b,t], -Fnl, pES[e,t,1])
    end
end
EV_2_sys(EVnode, Ratea, F, T, gpe, pfe, pEV) = for (a,node)=enumerate(EVnode), t=0:T
    JuMP.add_to_expression!(gpe[t], -1., pEV[a,t])
    for b=eachindex(Ratea)
        Fnl = F[b, node]
        JuMP.add_to_expression!(pfe[b,t], -Fnl, pEV[a,t])
    end
end
function add_pf_limits(m, T, CaD, Wnode, Ratea, pfe)
    JuMP.@variable(m, 0 <= pϵ[b=eachindex(Ratea), t=0:T])
    for (b,Lim)=enumerate(Ratea)
        d = CaD["branch"][string(b)]
        (d["f_bus"] ∈ Wnode || d["t_bus"] ∈ Wnode) && continue
        JuMP.@constraint(m, [t=0:T],  pfe[b,t] <= Lim + pϵ[b,t])
        JuMP.@constraint(m, [t=0:T], -Lim -pϵ[b,t] <= pfe[b,t])
    end
    pϵ
end

_η(rng) = rand(rng, 0.93:1e-4:0.97)
_1(ch, s, envs, a...) = (gs = envs[s]; Threads.@spawn(_3(ch, s, gs, a...)))
_3(ch, s, gs, a...) = (MIP_sub!(s, gs, a...); put!(ch, s))
function MIP_sub!(
        s, genv,
        Sb, S, T, Δtˈ1h_ratio,
        F, CaD, Ratea,
        Ggvec, Gnode, GPmax, GPref,
        LmuTup, Ltype, LPmax, Lnode,
        Wnode, WPmax, WmuMat, Wscnvec,
        EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
        EEmax, EEini, EPminDiv, Enode
    )
    rng = Random.Xoshiro(s)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variables(m, begin # common variables (binary + continuous)
        0 <= bES0[e=eachindex(EEmax), u=(0,1)] <= 1, Int
        0 <= bx0[a=eachindex(EVLmax)] <= 1, Int
        0 <= eV0[a=eachindex(EVEmax)] <= EVEmax[a]
        0 <= eS0[e=eachindex(EEmax)] <= EEmax[e]
    end)
    xs = [bES0...; bx0; eV0; eS0]
    N = length(bES0)+length(bx0)+length(eV0)+length(eS0)
    Pi=fill(0.,N); Xl=similar(Pi); Xl2=similar(Xl)
    # complete the rest
    bES = Dict((e,t,u) => t==0 ? bES0[e,u] : JuMP.@variable(m, lower_bound=0, upper_bound=1, integer=true) for e=eachindex(EEmax) for t=0:T for u=(0,1)) # here t really starts from 0
    bx = Dict((a,t) => t==0 ? bx0[a] : JuMP.@variable(m, lower_bound=0, upper_bound=1, integer=true) for a=eachindex(EVLmax) for t=-1:T)
    eV = Dict((a,t) => t==0 ? eV0[a] : JuMP.@variable(m, lower_bound=0, upper_bound=(1-max(t,0)/10T)EVEmax[a]) for a=eachindex(EVEmax) for t=-1:T)
    eS = Dict((e,t) => t==0 ? eS0[e] : JuMP.@variable(m, lower_bound=0, upper_bound=EEmax[e]) for e=eachindex(EEmax) for t=-1:T)
    # inner other variables
    JuMP.@variables(m, begin # Station variables
        0 <= bu[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= bv[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= pDR[a=eachindex(EVLmax), t=0:T, u=(0,1)] # DR of EVs (u=0: EV reduce their demand, u=1: EV increase their demand)
        pEV[a=eachindex(EVLmax), t=0:T]
    end)
    for (a,Pl)=enumerate(EVLmax) # Station inner
        JuMP.set_lower_bound(eV[a,T], max(0.2 * EVEmax[a], 0.4 * EVEini[a])) # End-of-horizon level
        JuMP.fix(bx[a,-1], Bxini[a]; force=true) # initial `x` variable
        JuMP.fix(eV[a,-1], EVEini[a]; force=true)
        L01vec = EVL01[a]
        ηc = _η(rng)
        for t=0:T
            JuMP.@constraints(m, begin
                bu[a,t]-bv[a,t] == bx[a,t]-bx[a,t-1]
                bx[a,t] ≤ 1-bv[a,t] # min down time is 1
                sum(bu[a,i] for i=range(t;step=-1,length=4+a) if i≥0) ≤ bx[a,t] # min up time is UT
                0.1bx[a,t]Pl ≤ pEV[a,t] # minimum charging power (timeless)
                pEV[a,t] ≤ 2.5bx[a,t]Pl # maximum charging power (timeless)
                eV[a,t]-eV[a,t-1] == ηc*(pEV[a,t]-L01vec[t+1]Pl)Δtˈ1h_ratio - pDR[a,t,1] + pDR[a,t,0]
            end)
        end
    end
    Qy = JuMP.AffExpr()
    JuMP.@variable(m, pES[e=eachindex(EEmax), t=0:T, u=(0,1)]) # this generates degradation cost: sum(pES)C_deg
    C01 = 1.5:0.1:6.5
    for (e,Tp)=enumerate(EPminDiv) # ES unit - inner constrs
        c01 = rand(rng, C01); for t=0:T, u=(0,1)
            JuMP.add_to_expression!(Qy, c01, pES[e,t,u]) # degrade cost
        end
        JuMP.@constraint(m, [t=0:T], bES[e,t,0] + bES[e,t,1] <= true) # only 3 states
        JuMP.@constraint(m, [t=0:T, u=(0,1)], (EEmax[e]/Tp[u+1])bES[e,t,u] <= pES[e,t,u]) # Min_power
        JuMP.@constraint(m, [t=0:T, u=(0,1)], pES[e,t,u] <= 2.0EEmax[e]bES[e,t,u]) # Max_power
        JuMP.set_lower_bound(eS[e,T], max(0.2 * EEmax[e], 0.5 * EEini[e]))
        JuMP.fix(eS[e,-1], EEini[e]; force=true)
        ηc, ηd = _η(rng), _η(rng)
        JuMP.@constraint(m, [t=0:T], eS[e,t]-eS[e,t-1] == (pES[e,t,1]ηc-pES[e,t,0]/ηd)Δtˈ1h_ratio)
    end

    pfe = Dict((b,t) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T)
    gpe = Dict(t => JuMP.AffExpr() for t=0:T)
    ES_2_sys(Enode, Ratea, F, T, gpe, pfe, pES)
    EV_2_sys(EVnode, Ratea, F, T, gpe, pfe, pEV)
    Cdr = 100:.01:500; for p=pDR JuMP.add_to_expression!(Qy, rand(rng, Cdr), p) end

    Wmat = Wscnvec[s] # WmuMat
    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T])
    Cϖ = 150:0.01:250; for ϖ=ϖ JuMP.add_to_expression!(Qy, rand(rng, Cϖ), ϖ) end
    for (w,node)=enumerate(Wnode) # wind - to system
        Pmax = WPmax[w]
        for t=0:T
            Pw = Pmax * Wmat[1+t, w]
            JuMP.set_upper_bound(ϖ[w,t], Pw)
            JuMP.add_to_expression!(gpe[t], -1., ϖ[w,t])
            JuMP.add_to_expression!(gpe[t], Pw)
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], -Fnl, ϖ[w,t])
                JuMP.add_to_expression!(pfe[b,t],  Fnl * Pw)
            end
        end
    end
    JuMP.@variable(m, 0 <= ζ[l=eachindex(LPmax), t=0:T])
    Cζ = 1000:.1:1800; for ζ=ζ JuMP.add_to_expression!(Qy, rand(rng, Cζ), ζ) end
    for (l,node)=enumerate(Lnode) # load - to system
        LdMax, L01Curve = LPmax[l], LmuTup[Ltype[l]]
        for t=0:T
            Pl = LdMax * L01Curve[1+t]
            JuMP.set_upper_bound(ζ[l,t], Pl)
            JuMP.add_to_expression!(gpe[t], ζ[l,t])
            JuMP.add_to_expression!(gpe[t], -Pl)
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], Fnl, ζ[l,t])
                JuMP.add_to_expression!(pfe[b,t], -Fnl * Pl)
            end
        end
    end
    JuMP.@variable(m, 0 <= r[g=eachindex(GPmax), t=0:T, u=(0,1)])
    Cr = 7:.1:37; for r=r JuMP.add_to_expression!(Qy, rand(rng, Cr), r) end
    for (g,node)=enumerate(Gnode) # Generator Redispatch
        Pmax = GPmax[g]
        for t=0:T
            Pba = Pmax * GPref[t+1,g]
            JuMP.set_upper_bound(r[g,t,0], Pba)
            JuMP.set_upper_bound(r[g,t,1], Pmax - Pba)
            JuMP.add_to_expression!(gpe[t], Pba)
            JuMP.add_to_expression!(gpe[t], -1., r[g,t,0])
            JuMP.add_to_expression!(gpe[t], r[g,t,1])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], Fnl * Pba)
                JuMP.add_to_expression!(pfe[b,t], -Fnl, r[g,t,0])
                JuMP.add_to_expression!(pfe[b,t], Fnl, r[g,t,1])
            end
        end
    end

    JuMP.@constraint(m, [t=0:T], gpe[t] == 0)
    pϵ = add_pf_limits(m, T, CaD, Wnode, Ratea, pfe);
    Cϵ = 9999:.1:11000; for ϵ=pϵ JuMP.add_to_expression!(Qy, rand(rng, Cϵ), ϵ) end
    
    Qy /= S # Prob
    JuMP.@objective(m, Min, Qy)

    Cd = Vector{Cdouble}(undef, N+2) # refers to indicies of _master_
    Sb[s] = (; m, o, refd, refi, xs, Qy, N, Xl, Pi, Xl2, Cd)
end

end

