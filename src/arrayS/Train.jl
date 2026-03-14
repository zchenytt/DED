module Train

import ..Settings, ..Build
import JuMP, Gurobi

function add_cut(n, #=Both for Fs&&mC: Obj=# Cn, xtrial, mst, mstLk, Ncuts,
    #=0.0 for a Fs, -1.0 for a mC=# θs_coeff::Cdouble, evt)
    o, om = n.moi_backend, mst.moi_backend
    ci, cd, st = n[:ci], n[:cd], n[:xlinkstart]
    len = Cint(length(xtrial))
    Gurobi.GRBgetdblattrarray(o, "RC", st, len, cd)
    cd[end] = θs_coeff
    Cn -= view(cd, 1:len)'xtrial
    cn = -Cn
    nnz = len + Cint(1)
    @lock(mstLk, Settings.addle(om, nnz, ci, cd, cn))
    Threads.atomic_add!(Ncuts, 1)
    notify(evt)
    nothing
end

function _7(s, tks, Ths, a...)
    c, S = -2, length(tks)
    while count(!istaskdone, tks) < Ths # `inn-tasks underoccupy hardware resources`
        t = tks[s]
        if istaskdone(t)
            wait(t)
            tks[s] = spawn_scene(s, a...) # must use function barrier here
        end
        s = mod(s, S)+1
        c === Ths && break
        c += 1
    end
    s
end
function _8(otk, a...)
    if istaskdone(otk)
        wait(otk)
        otk = Threads.@spawn(mst_procedure(a...)) # the bindings are const, thus save a function barrier
    end
    otk
end
wait3(stk, otr, tks) = (wait(stk); wait(otr.x); foreach(wait, tks))
ss(a...) = Threads.@spawn(schedulingloop(a...))
function schedulingloop(
    #=non_locals=# evt, s, otr,
    #=bookkeeping=# tks, Lk, ref, mst, inn, Ncuts,
    #=Settings=# THREAD_FOR_BLOCKS, MaxSec)
    otk, D_t, _ = otr.x, 1_000_000_000*MaxSec, notify(evt)
    t0 = time_ns()
    while true
        wait(evt) # irredundant as it sets status `false`
        d_t = time_ns() - t0
        # t_elapsed = d_t / 1_000_000_000
        # @ccall(printf("t_elapsed=%e\n"::Cstring; t_elapsed::Cdouble)::Cint)
        d_t > D_t && break
        otk = _8(otk, mst, Ncuts, Lk, ref, evt)
        s = _7(s, tks, THREAD_FOR_BLOCKS, inn, mst, Ncuts, Lk, ref, evt)
        if !evt.set
            sleep(0.001) # avoid a spinning loop
            notify(evt)
        end
    end
    otr.x = otk
end

function _3(m)
    Settings.opt_and_ter(m) === JuMP.OPTIMAL || error()
    Build.get_nt(m)
end
_2(k, lb) = @ccall(printf("k=%d, lb=%e\n"::Cstring; k::Cint, lb::Cdouble)::Cint)
function mst_procedure(m, Ncuts, Lk, ref, evt)
    Ncuts_out = Ncuts.value
    if m[:Ncuts] !== Ncuts_out
        nt = @lock(Lk.mst, _3(m))
        _2(nt.k, nt.lb) # Keep a Logging
        @lock(Lk.ref, setfield!(ref, :x, nt))
        notify(evt)
        m[:Ncuts] = Ncuts_out
    end
end

might_be_infeasible(t) = t === JuMP.INFEASIBLE || t === JuMP.INFEASIBLE_OR_UNBOUNDED
print_chgmod(s, str) = (@ccall(printf("s=%d >> %s\n"::Cstring; s::Cint, str::Cstring)::Cint); nothing)
print_vio(s, str, vio) = (@ccall(printf("s=%d-%s vio=%e\n"::Cstring; s::Cint, str::Cstring, vio::Cdouble)::Cint); nothing)

function analyze(lb, ub)
    agap = ub - lb
    rgap = agap / ub
    println("lb=$lb, agap=$agap, rgap=$rgap")
end
function post_check(tks, inn, #=mst=# m)
    Settings.opt_and_ter(m) === JuMP.OPTIMAL || error()
    nt = Build.get_nt(m) # This nt is not derived from `ref`
    ub, xtrial = Threads.Atomic{Float64}(nt.cˈx), nt.xlink
    for s = eachindex(tks)
        tks[s] = Threads.@spawn(write_to_ub(ub, inn, s, xtrial))
    end
    foreach(wait, tks)
    analyze(nt.lb, ub.value)
end
function write_to_ub(ub, inn, s, xtrial) # [post training] once for each block
    n, objval = inn[s], Inf
    if n[:ismC]
        Build.fix(n, xtrial)
        ter = Settings.opt_and_ter(n)
        if ter === JuMP.OPTIMAL
            objval = Settings.getmodeldblattr(n, "ObjVal")
        end
    end
    Threads.atomic_add!(ub, objval/length(inn))
end

spawn_scene(a...) = Threads.@spawn(scene_subprocedure(a...))
function scene_subprocedure(s, inn, #=mst=# m, Ncuts, Lk, ref, evt)
    n = inn[s]
    ismC = n[:ismC]
    nt = @lock(Lk.ref, ref.x) # This is static locally
    k_out, θs_check, xtrial = nt.k, nt.θ[s], nt.xlink
    if n[:k] !== k_out
        Build.fix(n, xtrial)
        n[:k] = k_out
    end
    ter = Settings.opt_and_ter(n)
    if ter === JuMP.OPTIMAL
        obj = Settings.getmodeldblattr(n, "ObjVal")
        if ismC
            vio = obj - θs_check
            # print_vio(s, "mC", vio)
            vio > 1e-5 && add_cut(n, obj, xtrial, m, Lk.mst, Ncuts, -1., evt)
            return nothing
        else # solved a Feas-OBJed-model to optimality
            vio = obj
            # print_vio(s, "Fs", vio)
            if vio > 1e-5
                add_cut(n, obj, xtrial, m, Lk.mst, Ncuts, 0., evt)
            else
                Build.chgabnd(n, 0) # chg mode to minC
            end
        end
    elseif ismC && might_be_infeasible(ter)
        Build.chgabnd(n, 1)
        # print_chgmod(s, "Fs")
    else
        error(20938475)
    end
end

end