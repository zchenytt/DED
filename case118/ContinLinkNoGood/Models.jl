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
const SB_CUT_COT = 1e-5
const SubLPTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refi::Base.RefValue{Int32}, refd::Base.RefValue{Float64}, Ci::Vector{Int32}, Cd::Vector{Float64}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, N::Int64}
const SubMIPTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refi::Base.RefValue{Int32}, refd::Base.RefValue{Float64}, Ci::Vector{Int32}, Cd::Vector{Float64}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, Bv::Vector{Int8}, N::Int64}

T2Δt(T) = if T === 16 || T === 12
        Δt_1h_ratio = 15/60
    elseif T === 48
        Δt_1h_ratio = 5/60
    else
        error("T = $T invalid")
end

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
function multi_θ_mst!(genv, S, Agi, EEmax)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variable(m, θ[1:S])                  # Notice the strict order
    JuMP.@variable(m, 0 ≤ p0[g=eachindex(Agi)] ≤ 1)    # [t=0]
    JuMP.@variable(m, 0 ≤ en0[e=eachindex(EEmax)] ≤ 1) # [t=0]
    N = length(p0) + length(en0)
    Θ, Xl, Xl2 = fill(0., S), fill(0., N), fill(0., N) # Notice the strict order
    printstyled("multi_θ master built with N = $N\n"; color = 30)
    Vv = CircularBuffer(fill(SB_CUT_COT, S))
    Σt = Ref(0.) # ReOptΣt
    (; m, o, refd, refi, S, N, Xl, Θ, Vv, Σt, Xl2)
end
function single_θ_mst!(genv, S, Agi, Agp, GPmax, KRmp, EEmax, EEini, EPmax, T)
    Δt_1h_ratio = T2Δt(T)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variable(m, θ)
    JuMP.@variable(m, p0[g=eachindex(Agi)])    # [t=0]
    JuMP.@variable(m, en0[e=eachindex(EEmax)]) # [t=0]
    for (e,Emax)=enumerate(EEmax)
        Em1, Pma, Emin, En0 = EEini[e], EPmax[e], 0., en0[e]
        JuMP.set_lower_bound(En0, max(Em1 - (Pma/0.95)Δt_1h_ratio, Emin))
        JuMP.set_upper_bound(En0, min(Em1 + (0.95*Pma)Δt_1h_ratio, Emax))
    end
    for (g,i)=enumerate(Agi)
        Pmin, Pmax, Pm1, Pg0 = 0., GPmax[i], Agp[g], p0[g]
        RmpLim = KRmp * Pmax
        JuMP.set_lower_bound(Pg0, max(Pm1 - RmpLim, Pmin))
        JuMP.set_upper_bound(Pg0, min(Pm1 + RmpLim, Pmax))
    end
    N = length(p0) + length(en0)
    Ci, Cd = Cint[1:N; 0], fill(0., N+2)
    Xl = fill(0., N) # Notice the strict order
    printstyled("single_θ master built with N = $N, T = $T\n"; color = 30)
    (; m, o, refd, refi, S, N, Xl, Ci, Cd)
end

"""
For the 2nd-stage cost, we additionally add a `c_s'z` bias term,
But the structure is almost unaltered.
The cut is pi*'x - theta_s < Gn*, where Gn* = pi*'x_check - Obj* (LP subproblem)
We still retrieve pi* via RC(z_copy), rather than introduce real `==`-constraints
"""
function LP_sub!(
        Ve, s, S, N, T, F, genv, Agi, Ratea,
        Lnode, Ltype, LmuTup, LPmax, LCshd,
        WPmax, Wscnvec, Wnode, WCcur,
        EEmax, EPmax, EEini, Enode, EClin, # EPmin is not considered in the LP model
        KRmp, GPmax, Gnode, GClin, KCres,
    )
    Δt_1h_ratio = T2Δt(T)
    m, o, ge, refi, refd = _9(genv)
    Ci, Cd = Cint[range(S; length=N); s-1], Vector{Cdouble}(undef, N+2) # This refers to indicies of _master_

    JuMP.@variable(m, p0[g=eachindex(Agi)])    # [t=0]
    JuMP.@variable(m, en0[e=eachindex(EEmax)]) # [t=0]
    Xl, Xl2, x_che = fill(0.,N), fill(0.,N), fill(0.,N) # x_che is locally static

    QyCz = JuMP.AffExpr()
    pfe = Dict((b,t) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T)
    gpe = Dict(t     => JuMP.AffExpr()                        for t=0:T)
    JuMP.@variables(m, begin # inner variables
        0 <= ϱ[g=eachindex(Agi), t=0:T] # effective generator power that injects into the network
        0 <= ζ[l=eachindex(Lnode), t=0:T]
        0 <= ϖ[w=eachindex(WPmax), t=0:T]
        0 <= en[e=eachindex(EEmax), t=1:T] <= EEmax[e]
        0 <= pe[e=eachindex(EPmax), t=0:T, u=(0,1)] <= EPmax[e]
        0 <= es0slk[e=eachindex(EPmax), u=(0,1)]
        0 <=  p[g=eachindex(Agi), t=1:T]
        0 <=  r[g=eachindex(Agi), t=0:T, u=(0,1)]
        -Ratea[b] <= pf[b=eachindex(Ratea), t=0:T] <= Ratea[b]
    end)

    # TODO Since we're currently using classic Benders cut, we don't tighten the region of the linking variable
    for (g,i)=enumerate(Agi) # generator - inner dynamics
        Pmax = GPmax[i]; PRmp = KRmp * Pmax
        JuMP.@constraint(m, p0[g] - r[g,0,0] >= 0)
        JuMP.@constraint(m, p0[g] + r[g,0,1] <= Pmax)
        JuMP.@constraint(m, [t=1:T], p[g,t] - r[g,t,0] >= 0)
        JuMP.@constraint(m, [t=1:T], p[g,t] + r[g,t,1] <= Pmax)
        JuMP.@constraint(m, p[g,1] <= p0[g] + PRmp)
        JuMP.@constraint(m, p[g,1] >= p0[g] - PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t] <= p[g,t-1] + PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t] >= p[g,t-1] - PRmp)
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,1]) >= 4)
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,0]) >= 2)
    end
    for e=eachindex(EPmax) # ES unit - inner dynamics
        t=0; JuMP.@constraint(m, EEini[e] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en0[e] + es0slk[e,1] - es0slk[e,0])
        t=1; JuMP.@constraint(m, en0[e]   + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
        JuMP.@constraint(m, [t=2:T], en[e,t-1] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
        JuMP.set_lower_bound(en[e,T], EEini[e])
    end
    for (e,node)=enumerate(Enode), t=0:T # ES unit - to system
        JuMP.add_to_expression!(QyCz, EClin, pe[e,t,0])
        JuMP.add_to_expression!(QyCz, EClin, pe[e,t,1])
        JuMP.add_to_expression!(gpe[t], pe[e,t,0])
        JuMP.add_to_expression!(gpe[t], -1., pe[e,t,1])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t],  Fnl, pe[e,t,0])
            JuMP.add_to_expression!(pfe[b,t], -Fnl, pe[e,t,1])
        end
    end
    for (g,i)=enumerate(Agi) # generator - to system
        gC, node = GClin[i], Gnode[i]
        for t=0:T JuMP.add_to_expression!(QyCz, KCres * gC, r[g,t,0]) end
        for t=0:T JuMP.add_to_expression!(QyCz, KCres * gC, r[g,t,1]) end
        JuMP.add_to_expression!(QyCz, gC, p0[g]); for t=1:T JuMP.add_to_expression!(QyCz, gC, p[g,t]) end
        t=0; JuMP.@constraint(m, ϱ[g,t] <= p0[g]); JuMP.add_to_expression!(gpe[t], ϱ[g,t])
        for t=1:T
            JuMP.@constraint(m, ϱ[g,t] <= p[g,t])
            JuMP.add_to_expression!(gpe[t], ϱ[g,t])
        end
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            t=0; JuMP.add_to_expression!(pfe[b,t], Fnl, ϱ[g,t])
            for t=1:T JuMP.add_to_expression!(pfe[b,t], Fnl, ϱ[g,t]) end
        end
    end
    for (l,node)=enumerate(Lnode) # load - to system
        ty = Ltype[l]
        for t=0:T
            Pmax, l01 = LPmax[l], LmuTup[ty][1+t]
            Pl = Pmax * l01
            JuMP.set_upper_bound(ζ[l,t], Pl)
            JuMP.add_to_expression!(gpe[t], -Pl)
            JuMP.add_to_expression!(gpe[t], ζ[l,t])
            JuMP.add_to_expression!(QyCz, LCshd, ζ[l,t])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], -Fnl * Pl)
                JuMP.add_to_expression!(pfe[b,t], Fnl, ζ[l,t])
            end
        end
    end
    Wmat = Wscnvec[s]; for (w,node)=enumerate(Wnode) # wind - to system
        Pmax = WPmax[w]
        for t=0:T
            w01 = Wmat[1+t, w]
            Pw = Pmax * w01
            JuMP.set_upper_bound(ϖ[w,t], Pw)
            JuMP.add_to_expression!(gpe[t], Pw)
            JuMP.add_to_expression!(gpe[t], -1., ϖ[w,t])
            JuMP.add_to_expression!(QyCz, WCcur, ϖ[w,t])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t],  Fnl * Pw)
                JuMP.add_to_expression!(pfe[b,t], -Fnl, ϖ[w,t])
            end
        end
    end
    for v=es0slk JuMP.add_to_expression!(QyCz, 2000., v) end
    JuMP.@objective(m, Min, QyCz); JuMP.@constraints(m, begin
        [b=eachindex(Ratea), t=0:T], pfe[b,t] == pf[b,t]
        [t=0:T], gpe[t] == 0
    end)

    Ve[s] = (; m, o, refi, refd, Ci, Cd, Xl, Xl2, x_che, N)
end

function MIP_sub!(
        Ve, s, S, N, T, F, genv, Agi, Agp, Ratea,
        Lnode, Ltype, LmuTup, LPmax, LCshd,
        WPmax, Wscnvec, Wnode, WCcur,
        EEmax, EPmax, EEini, Enode, EPmin,
        KRmp, GPmax, Gnode, GClin, KCres,
    )
    Δt_1h_ratio = T2Δt(T)
    m, o, ge, refi, refd = _9(genv)
    Ci, Cd = Cint[range(S; length=N); s-1], Vector{Cdouble}(undef, N+2) # This refers to indicies of _master_
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "Presolve", 2)
    # (All Continuous) Linking Normalized Variables
    JuMP.@variable(m, p0N[g=eachindex(Agi)])    # [t=0]
    JuMP.@variable(m, en0N[e=eachindex(EEmax)]) # [t=0]
    #=Pi=# Xl,#=Aux=# Xl2,#=LocalStatic=# x_che = fill(0.,N),fill(0.,N),fill(0.,N)
    #= Inner Binaries =# JuMP.@variable(m, 0 <= be[e=eachindex(EPmax), t=0:T, u=(0,1)] <= 1)
    Bv = fill(Cchar('C'), length(be))

    JuMP.@variable(m, p0[g=eachindex(Agi)])   
    JuMP.@variable(m, en0[e=eachindex(EEmax)])
    for (g,i)=enumerate(Agi)
        Pmin, Pmax, Pm1 = 0., GPmax[i], Agp[g]
        RmpLim = KRmp * Pmax
        lb = max(Pm1 - RmpLim, Pmin)
        gap = min(Pm1 + RmpLim, Pmax) - lb 
        JuMP.@constraint(m, p0[g] == lb + gap * p0N[g])
    end
    for (e,Emax)=enumerate(EEmax)
        Em1, Pma, Emin  = EEini[e], EPmax[e], 0.
        lb = max(Em1 - (Pma/0.95)Δt_1h_ratio, Emin)
        gap = min(Em1 + (0.95*Pma)Δt_1h_ratio, Emax) - lb
        JuMP.@constraint(m, en0[e] == lb + gap * en0N[e])
    end

    #= Inner Else =# Qy = JuMP.AffExpr()
    pfe = Dict((b,t) => JuMP.AffExpr() for b=eachindex(Ratea) for t=0:T)
    gpe = Dict(t     => JuMP.AffExpr()                        for t=0:T)
    JuMP.@variables(m, begin # inner variables
        0 <= ϱ[g=eachindex(Agi), t=0:T] # effective generator power that injects into the network
        0 <= ζ[l=eachindex(Lnode), t=0:T]
        0 <= ϖ[w=eachindex(WPmax), t=0:T]
        0 <= en[e=eachindex(EEmax), t=1:T] <= EEmax[e]
        pe[e=eachindex(EPmax), t=0:T, u=(0,1)]
        0 <= es0slk[e=eachindex(EPmax), u=(0,1)]
        0 <=  p[g=eachindex(Agi), t=1:T]
        0 <=  r[g=eachindex(Agi), t=0:T, u=(0,1)]
        -Ratea[b] <= pf[b=eachindex(Ratea), t=0:T] <= Ratea[b]
    end)
    for (g,i)=enumerate(Agi) # generator - inner dynamics
        Pmax = GPmax[i]; PRmp = KRmp * Pmax
        JuMP.@constraint(m, p0[g] - r[g,0,0] >= 0)
        JuMP.@constraint(m, p0[g] + r[g,0,1] <= Pmax)
        JuMP.@constraint(m, [t=1:T], p[g,t] - r[g,t,0] >= 0)
        JuMP.@constraint(m, [t=1:T], p[g,t] + r[g,t,1] <= Pmax)
        JuMP.@constraint(m, p[g,1] <= p0[g] + PRmp)
        JuMP.@constraint(m, p[g,1] >= p0[g] - PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t] <= p[g,t-1] + PRmp)
        JuMP.@constraint(m, [t=2:T], p[g,t] >= p[g,t-1] - PRmp)
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,1]) >= 4)
        JuMP.@constraint(m, [t=T÷2:T], sum(r[:,t,0]) >= 2)
    end
    for e=eachindex(EPmax) # ES unit - inner dynamics
        JuMP.@constraint(m, [t=0:T], be[e,t,0] + be[e,t,1] <= 1)
        JuMP.@constraint(m, [t=0:T, u=(0,1)], be[e,t,u]EPmin[e] <= pe[e,t,u])
        JuMP.@constraint(m, [t=0:T, u=(0,1)], pe[e,t,u] <= be[e,t,u]EPmax[e])
        t=0; JuMP.@constraint(m, EEini[e] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en0[e] + es0slk[e,1] - es0slk[e,0])
        t=1; JuMP.@constraint(m, en0[e]   + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
        JuMP.@constraint(m, [t=2:T], en[e,t-1] + (0.95pe[e,t,1] - pe[e,t,0]/0.95)Δt_1h_ratio == en[e,t])
        JuMP.set_lower_bound(en[e,T], EEini[e])
    end
    for (e,node)=enumerate(Enode), t=0:T # ES unit - to system
        JuMP.add_to_expression!(gpe[t], pe[e,t,0])
        JuMP.add_to_expression!(gpe[t], -1., pe[e,t,1])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t],  Fnl, pe[e,t,0])
            JuMP.add_to_expression!(pfe[b,t], -Fnl, pe[e,t,1])
        end
    end
    for (g,i)=enumerate(Agi) # generator - to system
        gC, node = GClin[i], Gnode[i]
        for t=0:T JuMP.add_to_expression!(Qy, KCres * gC, r[g,t,0]) end
        for t=0:T JuMP.add_to_expression!(Qy, KCres * gC, r[g,t,1]) end
        JuMP.add_to_expression!(Qy, gC, p0[g]); for t=1:T JuMP.add_to_expression!(Qy, gC, p[g,t]) end
        t=0; JuMP.@constraint(m, ϱ[g,t] <= p0[g]); JuMP.add_to_expression!(gpe[t], ϱ[g,t])
        for t=1:T
            JuMP.@constraint(m, ϱ[g,t] <= p[g,t])
            JuMP.add_to_expression!(gpe[t], ϱ[g,t])
        end
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            t=0; JuMP.add_to_expression!(pfe[b,t], Fnl, ϱ[g,t])
            for t=1:T JuMP.add_to_expression!(pfe[b,t], Fnl, ϱ[g,t]) end
        end
    end
    for (l,node)=enumerate(Lnode) # load - to system
        ty = Ltype[l]
        for t=0:T
            Pmax, l01 = LPmax[l], LmuTup[ty][1+t]
            Pl = Pmax * l01
            JuMP.set_upper_bound(ζ[l,t], Pl)
            JuMP.add_to_expression!(gpe[t], -Pl)
            JuMP.add_to_expression!(gpe[t], ζ[l,t])
            JuMP.add_to_expression!(Qy, LCshd, ζ[l,t])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t], -Fnl * Pl)
                JuMP.add_to_expression!(pfe[b,t], Fnl, ζ[l,t])
            end
        end
    end
    Wmat = Wscnvec[s]; for (w,node)=enumerate(Wnode) # wind - to system
        Pmax = WPmax[w]
        for t=0:T
            w01 = Wmat[1+t, w]
            Pw = Pmax * w01
            JuMP.set_upper_bound(ϖ[w,t], Pw)
            JuMP.add_to_expression!(gpe[t], Pw)
            JuMP.add_to_expression!(gpe[t], -1., ϖ[w,t])
            JuMP.add_to_expression!(Qy, WCcur, ϖ[w,t])
            for b=eachindex(Ratea)
                Fnl = F[b, node]
                JuMP.add_to_expression!(pfe[b,t],  Fnl * Pw)
                JuMP.add_to_expression!(pfe[b,t], -Fnl, ϖ[w,t])
            end
        end
    end
    for v=es0slk JuMP.add_to_expression!(Qy, 2000., v) end
    JuMP.@objective(m, Min, Qy); JuMP.@constraints(m, begin
        [b=eachindex(Ratea), t=0:T], pfe[b,t] == pf[b,t]
        [t=0:T], gpe[t] == 0
    end)

    Ve[s] = (; m, o, refi, refd, Ci, Cd, Xl, Xl2, x_che, Bv, N)
end

end

