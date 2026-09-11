module ADMM

const W_fg, W_bg = Threads.nthreads(:interactive), Threads.nthreads(:default)
const ρ = 9e3
import ..Settings, JuMP, Gurobi, Random

function post_adjust_λi!(λi)
    r, S = sum(λi), length(λi)
    r ./= S
    for (s, λ) = enumerate(λi)
        @. λi[s] = λ - r
    end
end

upd_xi_objs!(admv, λi, xi_x) = for (s, λi) = enumerate(λi)
    m, X = admv[s], xi_x[end]
    fi, xi = m[:psQy], m[:xi] # decision vector
    JuMP.@objective(m, Min, fi + λi'xi + ρ/2 * sum((a-b)^2 for (a,b)=zip(xi,X)))
end
post_upd_xi_objs!(admv, λi, xi_x) = for (s, λi) = enumerate(λi)
    m = admv[s]
    fi, xi = m[:psQy], m[:xi] # decision vector
    JuMP.@objective(m, Min, fi + λi'xi)
end

update_dual!(λi, xi_x) = for (s, λold) = enumerate(λi)
    @. λi[s] = λold + ρ * (xi_x[s]-xi_x[end])
end

function upd_x!(xi_x, λi)
    m = Settings.Model(Settings.Env())
    l = length(xi_x[end])
    JuMP.@variable(m, x[1:l])
    JuMP.@objective(m, Min, (ADMM.ρ/2)sum(sum((a-b)^2 for (a,b)=zip(x,y)) for y=xi_x[1:end-1]) - x'sum(λi))
    JuMP.optimize!(m)
    JuMP.termination_status(m) == JuMP.OPTIMAL || error()
    @. xi_x[end] = JuMP.value(x)
end

function solve_all_and_check(admv)
    Threads.@threads for s=eachindex(admv)
        JuMP.optimize!(admv[s])
    end
    for s=eachindex(admv)
        ter = JuMP.termination_status(admv[s])
        (ter == JuMP.OPTIMAL || ter == JuMP.TIME_LIMIT) || error()
    end
end

upd_all_xi!(xi_x, admv) = for s=eachindex(admv)
    x = admv[s][:xi]
    @. xi_x[s] = JuMP.value(x)
end

function get_residual(xi_x)
    x = xi_x[end]
    maximum(maximum(abs(a-b) for (a,b)=zip(x,y)) for y=xi_x[1:end-1])
end

function _9(genv)
    m = Settings.Model(genv)
    o, refi, refd = m.moi_backend, Ref{Cint}(), Ref{Cdouble}()
    ge = Gurobi.GRBgetenv(o)
    m, o, ge, refi, refd
end
_η(rng) = rand(rng, 0.93:1e-4:0.97)
function Model(
    genv, s, S, T, Δtˈ1h_ratio,
    Ratea, F, CaD,
    EVnode, EVLmax, EVL01, EVEmax, EVEini,
    Bxini, Enode, EEmax, EEini, EPminDiv,
    WPmax,Wscnvec,Wnode,
    LPmax,LmuTup,Ltype,Lnode,
    GPmax,GPref,Gnode
)
    Qy, rng = JuMP.AffExpr(), Random.Xoshiro(s)
    m, o, ge, refi, refd = _9(genv)
    JuMP.@variables(m, begin # a superset of the complicating variables
        bES[e=eachindex(EEmax), t=0:T, u=(0,1)], Bin
        bx[a=eachindex(EVLmax), t=-1:T], Bin
        0 <= eS[e=eachindex(EEmax), t=-1:T] <= EEmax[e]
        0 <= eV[a=eachindex(EVEmax), t=-1:T] <= (1-max(t,0)/10T)EVEmax[a]
    end)
    m[:xi] = [bES[:,0,:]..., bx[:,0]..., eS[:,0]..., eV[:,0]...] # complicating variables

    JuMP.@variables(m, begin # Station variables
        0 <= bu[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= bv[a=eachindex(EVLmax), t=0:T] <= 1
        0 <= pDR[a=eachindex(EVLmax), t=0:T, u=(0,1)]
        pEV[a=eachindex(EVLmax), t=0:T]
    end)
    for (a,Pl)=enumerate(EVLmax) # Station inner
        JuMP.set_lower_bound(eV[a,T], max(.2 * EVEmax[a], .4 * EVEini[a])) # End-of-horizon level
        JuMP.unset_binary(bx[a,-1]); JuMP.fix(bx[a,-1], Bxini[a]; force=true) # initial `x` variable
        JuMP.fix(eV[a,-1], EVEini[a]; force=true)
        L01vec, ηc = EVL01[a], _η(rng)
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
    
    JuMP.@variable(m, pES[e=eachindex(EEmax), t=0:T, u=(0,1)])
    C01 = 1.5:0.1:6.5 # degradation cost
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
    for (e,node)=enumerate(Enode), t=0:T
        JuMP.add_to_expression!(gpe[t], pES[e,t,0])
        JuMP.add_to_expression!(gpe[t], -1., pES[e,t,1])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t],  Fnl, pES[e,t,0])
            JuMP.add_to_expression!(pfe[b,t], -Fnl, pES[e,t,1])
        end
    end
    for (a,node)=enumerate(EVnode), t=0:T
        JuMP.add_to_expression!(gpe[t], -1., pEV[a,t])
        for b=eachindex(Ratea)
            Fnl = F[b, node]
            JuMP.add_to_expression!(pfe[b,t], -Fnl, pEV[a,t])
        end
    end
    Cdr = 100:.01:500
    for p=pDR JuMP.add_to_expression!(Qy, rand(rng, Cdr), p) end

    JuMP.@variable(m, 0 <= ϖ[w=eachindex(WPmax), t=0:T])
    Cϖ = 150:0.01:250
    for p=ϖ JuMP.add_to_expression!(Qy, rand(rng, Cϖ), p) end
    for (w,node)=enumerate(Wnode) # wind - to system
        Pmax, Wmat = WPmax[w], Wscnvec[s] # WmuMat
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
    Cζ = 1000:.1:1800
    for p=ζ JuMP.add_to_expression!(Qy, rand(rng, Cζ), p) end
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
    Cr = 7:.1:37
    for p=r JuMP.add_to_expression!(Qy, rand(rng, Cr), p) end
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
    JuMP.@variable(m, 0 <= pϵ[b=eachindex(Ratea), t=0:T])
    for (b,Lim)=enumerate(Ratea)
        d = CaD["branch"][string(b)]
        (d["f_bus"] ∈ Wnode || d["t_bus"] ∈ Wnode) && continue
        JuMP.@constraint(m, [t=0:T],  pfe[b,t] <= Lim + pϵ[b,t])
        JuMP.@constraint(m, [t=0:T], -Lim -pϵ[b,t] <= pfe[b,t])
    end
    Cϵ = 9999:.1:11000
    for p=pϵ JuMP.add_to_expression!(Qy, rand(rng, Cϵ), p) end
    m[:psQy] = Qy/S
    JuMP.@objective(m, Min, m[:psQy])
    m
end

end
