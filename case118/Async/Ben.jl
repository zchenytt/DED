module Ben
import ..Models, ..Settings, Gurobi

const W_calc = Threads.nthreads()
const W_modeling = min((W_calc+1)÷2, W_calc)
const COT = 0.1

function load_che(mst)
    (; o, S, N, Θ, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 0, S,  Θ)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
end
function _fix_sub(n, x_che)
    (; N, o) = n
    Gurobi.GRBsetdblattrarray(o, "LB", 0, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 0, N, x_che)
end
function _load_RC(n)
    (; N, o, Xl) = n
    Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Xl)
end
function set_mst_θ_obj(mst)
    (; o, Θ, S) = mst
    # ✅ Note: since we had let `c'x === 0`, Minimize mean(Θ) == Minimize sum(Θ), the solution (x_che, Θs_che) should be the same
    Θ .= 1; Gurobi.GRBsetdblattrarray(o, "Obj", 0, S, Θ)
end

function eval_gap(mst, sub, S)
    Settings.opt_ass_opt(mst, "master")
    load_che(mst)
    lb = Settings.getmodeldblattr(mst, "ObjVal") / S
    x_che = mst.Xl
    Threads.@threads for n = sub
        _fix_sub(n, x_che)
    end
    Threads.@threads for n = sub
        Settings.opt_ass_opt(n)
    end
    ub = 0.
    for n = sub
        ObjVal = Settings.getmodeldblattr(n, "ObjVal")
        ub += ObjVal/S
    end
    agap = ub - lb
    rgap = agap / ub
    printstyled("lb = $lb, ub = $ub, rgap = $rgap, agap = $agap\n"; color = 21)
end

function _2(n, mst)
    (; o, N, Xl) = mst
    _fix_sub(n, Xl)
    Settings.opt_ass_opt(n)
    _load_RC(n)
    ObjVal = Settings.getmodeldblattr(n, "ObjVal")
    n.Cd[1:N] .= n.Xl
    n.Cd[N+1]  = -1
    n.Cd[end]  = Gn = n.Xl'Xl - ObjVal
end
_3(ch, s, sub, mst) = (_2(sub[s], mst); put!(ch, s))
function warm(mst, sub, S)
    if S > W_calc
        ch = Channel{Int}()
        for s=1:W_calc Threads.@spawn(_3(ch, s, sub, mst)) end
        i = W_calc
        while true
            i+1 > S && break
            take!(ch)
            s = i+=1
            Threads.@spawn(_3(ch, s, sub, mst))
        end
        for _=1:W_calc take!(ch) end
    else
        Threads.@threads for s=eachindex(sub)
            _2(sub[s], mst)
        end
    end
    for n=sub Gurobi.GRBaddconstr(mst.o, mst.N+1, n.Ci, n.Cd, Cchar('<'), n.Cd[end], C_NULL) end
end

function _4(ch, s, n)
    x_che, pi, N, Cd = n.x_che, n.Xl, n.N, n.Cd
    _fix_sub(n, x_che)
    Settings.opt_ass_opt(n, string("sub", s))
    _load_RC(n)
    Cd[1:N] .= pi
    Cd[N+1]  = -1
    Cd[end]  = Gn = pi'x_che - Settings.getmodeldblattr(n, "ObjVal")
    put!(ch, s)
end
_5(ch, CanSpawn, s, x_che, n) = (n.x_che .= x_che; CanSpawn[s]=false; Threads.@spawn(_4(ch, s, n)))
function main(mst, sub, S, Seconds)
    @assert S > W_calc
    Settings.opt_ass_opt(mst, "master")
    load_che(mst)
    x_che, Θ, N = mst.Xl, mst.Θ, mst.N
    ch, CanSpawn = Channel{Int}(), trues(S)
    for s=1:W_calc _5(ch, CanSpawn, s, x_che, sub[s]) end
    iun, ȷ, scn, ȷΣ, lb = W_calc, 0, 0, 0, -Inf
    t0 = time_ns()
    while true
        if !isready(ch) || scn > S
            if ȷ > 0
                Settings.opt_ass_opt(mst, "master") # update
                load_che(mst)
                lb = Settings.getmodeldblattr(mst, "ObjVal")
                ȷΣ += ȷ
                # @ccall(printf("ȷ=%d, ȷΣ=%d, lb=%.3e\n"::Cstring; ȷ::Cint, ȷΣ::Cint, lb::Cdouble)::Cint)
            else
                t_elapsed = 1e-9(time_ns() - t0)
                # @ccall(printf("t_elapsed=%.1e\n"::Cstring; t_elapsed::Cdouble)::Cint)
            end
            ȷ = scn = 0
            1e-9(time_ns() - t0) > Seconds && break
        end
        s = take!(ch)
        CanSpawn[s]=true; n=sub[s]; Θs_che=Θ[s]; Cd=n.Cd; Gn=Cd[end]
        vio = view(Cd, 1:N)'x_che - Gn - Θs_che
        if vio > COT
            # @ccall(printf("s = %d, vio = %.3e\n"::Cstring; s::Cint, vio::Cdouble)::Cint)
            Gurobi.GRBaddconstr(mst.o, N+1, n.Ci, Cd, Cchar('<'), Gn, C_NULL)
            ȷ += 1
        end
        while true
            iun += 1; s = Settings.to_1S(iun, S)
            CanSpawn[s] && (_5(ch, CanSpawn, s, x_che, sub[s]); break)
        end
        scn += 1
    end
end

"Build sub::Vector in parallel"
_1(ch, s, sub, S, a...) = (Models.sub!(sub, s, S, a...); put!(ch, s))
sub!(sub, S, a...) = if S > W_modeling
    ch = Channel{Int}()
    for s=1:W_modeling Threads.@spawn(_1(ch, s, sub, S, a...)) end
    i = W_modeling
    while true
        i+1 > S && break
        take!(ch)
        s = i+=1
        Threads.@spawn(_1(ch, s, sub, S, a...))
    end
    for _=1:W_modeling take!(ch) end
else
    Threads.@threads(for s=1:S Models.sub!(sub, s, S, a...) end)
end

end
