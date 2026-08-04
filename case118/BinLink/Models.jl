"""
sub: 1. LP; 2. MIP
mst: 1. single theta; 2. multi theta

The c'x term originally in 1st-stage can be put off to 2nd-stage.
This will have the benefit, if 1st-stage cost is nonlinear, then it
can be absorbed into the cut, such that our master problem is always LP.
"""
module Models
import DataStructures.CircularBuffer
import ..Settings, JuMP, Gurobi
const SB_CUT_COT = 0.03
const SubMIPTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, Bv::Vector{Int8}, N::Int64, Pi::Vector{Float64}, Xl::Vector{Float64}, Xl2::Vector{Float64}, Ci::Vector{Int32}, Cd::Vector{Float64}}

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
function multi_θ_mst!(genv, S, EVLmax, EEmax)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variable(m, θ[1:S]) # Notice the strict order
    Θ = fill(0.,S)
    JuMP.@variable(m, 0 <= bx0[a=eachindex(EVLmax)] <= 1)
    JuMP.@variable(m, 0 <= bES0[e=eachindex(EEmax), u=(0,1)] <= 1)
    N = length(bx0)+length(bES0); Xl=fill(0.,N); Bv=fill(Cchar('C'), N)
    JuMP.@constraint(m, [e=eachindex(EEmax)], bES0[e,0] + bES0[e,1] <= true)
    printstyled("multi_θ master built with N = $N\n"; color = 30)
    Σt = Ref(0.) # ReOptΣt
    Vv = CircularBuffer(fill(SB_CUT_COT, S))
    (; m, o, refd, refi, S, N, Xl, Θ, Bv, Σt, Vv)
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
_01v(m) = JuMP.@variable(m, lower_bound = 0, upper_bound = 1)
bx_scalar(m, a, t, bx0) = t==0 ? bx0[a] : _01v(m)
bES_scalar(m, e, t, u, bES0) = t==0 ? bES0[e,u] : _01v(m)

function MIP_sub!(
    Sb, s, S, genv, T, Δtˈ1h_ratio,
    F, CaD, Ratea,
    Wnode, WPmax, WmuMat, Wscnvec,
    LmuTup, Ltype, LPmax, Lnode,
    EEmax, EEini, EPminDiv, Enode,
    EVL01, EVLmax, EVEmax, EVEini, Bxini, EVnode,
    )
    m, o, ge, refi, refd = _9(genv)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "Presolve", 2)
    # linking variables (pure Bin)
    JuMP.@variable(m, 0 <= bx0[a=eachindex(EVLmax)] <= 1)
    JuMP.@variable(m, 0 <= bES0[e=eachindex(EEmax), u=(0,1)] <= 1)
    N = length(bx0)+length(bES0); Pi=fill(0.,N); Xl=similar(Pi); Xl2=similar(Xl); 
    # inner Bin variables
    bx = Dict((a,t) => bx_scalar(m, a, t, bx0) for a=eachindex(EVLmax) for t=-1:T)
    bES = Dict((e,t,u) => bES_scalar(m, e, t, u, bES0) for e=eachindex(EEmax) for t=0:T for u=(0,1)) # here t really starts from 0
    Bv = fill(Cchar('C'), length(bx)+length(bES)) # a long vector, the inner Bin length (N2B) can be deduced
    # inner other variables
    JuMP.@variables(m, begin # Station variables
        0 <= bu[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= bv[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= pDR[a=eachindex(EVLmax), t=0:T, u=(0,1)] # DR of EVs (u=0: EV reduce their demand, u=1: EV increase their demand)
        0 <= eV[a=eachindex(EVEmax), t=-1:T] <= (1-max(t,0)/10T)EVEmax[a] # storage UB (time variant)
        pEV[a=eachindex(EVLmax), t=0:T]
    end)
    for (a,Pl)=enumerate(EVLmax) # Station inner
        JuMP.set_lower_bound(eV[a,T], max(0.2 * EVEmax[a], 0.4 * EVEini[a])) # End-of-horizon level
        JuMP.fix(bx[a,-1], Bxini[a]; force=true) # initial `x` variable
        JuMP.fix(eV[a,-1], EVEini[a]; force=true)
        L01vec = EVL01[a]
        for t=0:T
            JuMP.@constraints(m, begin
                bu[a,t]-bv[a,t] == bx[a,t]-bx[a,t-1]
                bx[a,t] ≤ 1-bv[a,t] # min down time is 1
                sum(bu[a,i] for i=range(t;step=-1,length=4+a) if i≥0) ≤ bx[a,t] # min up time is UT
                0.1bx[a,t]Pl ≤ pEV[a,t] # minimum charging power (timeless)
                pEV[a,t] ≤ 2.5bx[a,t]Pl # maximum charging power (timeless)
                eV[a,t]-eV[a,t-1] == .95(pEV[a,t]-L01vec[t+1]Pl)Δtˈ1h_ratio - pDR[a,t,1] + pDR[a,t,0]
            end)
        end
    end

    JuMP.@variable(m, 0 <= eS[e=eachindex(EEmax), t=-1:T] <= EEmax[e])
    JuMP.@variable(m, pES[e=eachindex(EEmax), t=0:T, u=(0,1)]) # this generates degradation cost: sum(pES)C_deg
    for (e,Tp)=enumerate(EPminDiv) # ES unit - inner constrs
        JuMP.@constraint(m, [t=0:T], bES[e,t,0] + bES[e,t,1] <= true) # only 3 states
        JuMP.@constraint(m, [t=0:T, u=(0,1)], (EEmax[e]/Tp[u+1])bES[e,t,u] <= pES[e,t,u]) # Min_power
        JuMP.@constraint(m, [t=0:T, u=(0,1)], pES[e,t,u] <= 2.0EEmax[e]bES[e,t,u]) # Max_power
        JuMP.set_lower_bound(eS[e,T], max(0.2 * EEmax[e], 0.5 * EEini[e]))
        JuMP.fix(eS[e,-1], EEini[e]; force=true)
        JuMP.@constraint(m, [t=0:T], eS[e,t]-eS[e,t-1] == (.95pES[e,t,1]-pES[e,t,0]/.95)Δtˈ1h_ratio)
    end

    pfe = Dict((b,t) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T)
    gpe = Dict(t => JuMP.AffExpr() for t=0:T)
    ES_2_sys(Enode, Ratea, F, T, gpe, pfe, pES)
    EV_2_sys(EVnode, Ratea, F, T, gpe, pfe, pEV)

    Wmat = Wscnvec[s] # WmuMat
    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T])
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
    for (l,node)=enumerate(Lnode) # load - to system
        LdMax, L01Curve = LPmax[l], LmuTup[Ltype[l]]
        for t=0:T
            Pl = LdMax * L01Curve[1+t]
            JuMP.add_to_expression!(gpe[t], -Pl)
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], -Fnl * Pl)
            end
        end
    end
    JuMP.@constraint(m, [t=0:T], gpe[t] == 0)
    pϵ = add_pf_limits(m, T, CaD, Wnode, Ratea, pfe)
    JuMP.@objective(m, Min, sum(pES) + sum(pDR) + sum(ϖ) + 100sum(pϵ))
    
    Ci, Cd = Cint[range(S; length=N); s-1], Vector{Cdouble}(undef, N+2) # refers to indicies of _master_
    Sb[s] = (; m, o, refd, refi, Bv, N, Pi, Xl, Xl2, Ci, Cd) # Xl used to store local static x_che
end

end

# function EV_charging_station() 
#     JuMP.@variable(m, x[t=0:T], Bin) # 🟥 1: charging, 0: idle
#     JuMP.@variable(m, 0 ≤ u[t=0:T] ≤ 1) # jump to x=1
#     JuMP.@variable(m, 0 ≤ v[t=0:T] ≤ 1) # fall to x=0
#     JuMP.@constraint(m, x[0] - x_ini == u[0] - v[0]) # logic constr
#     JuMP.@constraint(m, [t=1:T], x[t]-x[t-1] == u[t]-v[t]) # logic constr
#     JuMP.@constraint(m, [t=0:T], x[t] ≤ 1 - v[t]) # min down time is 1
#     JuMP.@constraint(m, [t=0:T], sum(u[i] for i=range(t; step=-1, length=UT) if i ≥ 0) ≤ x[t]) # min up time is UT
    
#     JuMP.@variable(m, pEVc[t=0:T]) # charging power (defined in terms of the network side)
#     JuMP.@constraint(m, [t=0:T], x[t]pEVcmin ≤ pEVc[t])
#     JuMP.@constraint(m, [t=0:T], pEVc[t] ≤ x[t]pEVcmax)

#     JuMP.@variable(m, 0 ≤ eEV[t=0:T] ≤ eEVmax[t]) # 🟦 (time-dependent) electricity level
#     JuMP.@constraint(m, eEV[0]-eEV_ini == (Δt_1h_ratio)pEVc[0]ηEV - LoadEV[0]) # ηEV = 0.95
#     JuMP.@constraint(m, [t=1:T], eEV[t]-eEV[t-1] == (Δt_1h_ratio)pEVc[t]ηEV - LoadEV[t])
#     JuMP.set_lower_bound(eEV[T], eEV_ini) # this avoid emptying the electricity at end-of-horizon 
#     # The EV charging station's effect on the IEEE118 network is via `pEVc` (deemed a flexible load)
# end
# function ES_unit() # e.g. CAES
#     JuMP.@variable(m, 0 <= be[e=eachindex(EPmax), t=0:T, u=(0,1)] <= 1, Bin) # 🟥
#     JuMP.@constraint(m, [t=0:T], be[e,t,0] + be[e,t,1] <= 1) # 🟥

#     JuMP.@variable(m, pe[e=eachindex(EPmax), t=0:T, u=(0,1)])
#     JuMP.@constraint(m, [t=0:T, u=(0,1)], be[e,t,u]EPmin[e] <= pe[e,t,u])
#     JuMP.@constraint(m, [t=0:T, u=(0,1)], pe[e,t,u] <= be[e,t,u]EPmax[e])

#     JuMP.@variable(m, 0 <= en[e=eachindex(EEmax), t=1:T] <= EEmax[e]) # 🟦
#     t=0; JuMP.@constraint(m, EEini[e] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en0[e] + es0slk[e,1] - es0slk[e,0])
#     t=1; JuMP.@constraint(m, en0[e]   + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
#     JuMP.@constraint(m, [t=2:T], en[e,t-1] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
#     JuMP.set_lower_bound(en[e,T], EEini[e])
# end


