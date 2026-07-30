"""
Experiment result indicates that Classic Benders cut alone is not sufficient to close the Gap for MIP recourse subproblems (e.g. rgap can only be reduced to 19.5% and no more improvement)
"""
# TODO: Linking variable shouldn't be continuous if recourse is MIP---The primal ub would have a significant gap above the global lower bound. (while if linking variables are binary, primal ub wouldn't jump) 
# (Yes): See if (Ben-and then->SBen) cut can help reduce the gap

module MultiTheta
import ..Models, ..Settings, Gurobi
const ϵ = Models.SB_CUT_COT # not being used by single-theta classic Benders
const LP_CUT_OFF = 0.9 # If Classic Benders Cut can cut off more than This value, then add Classic Benders Cut directly
const Wcpu = Sys.CPU_THREADS - 1
const WcpuEval = Sys.CPU_THREADS ÷ 2
const EvalTimeLim = 60.
const EvalMIPGap = 0.
const TrainTimeLim = 24.
const TrainMIPGap = 9e-4

function _8(mst, Σt, ȷ, te)
    (; S, Xl, Xl2) = mst
    Settings.opt_ass_opt(mst, "master LP reOpt")
    load_che(mst)
    ReOptΣt = Σt.x += tReopt = Settings.getmodeldblattr(mst, "Runtime")
    ȷμ = Settings.getmodelintattr(mst, "NumConstrs")/S
    @ccall(printf("▶ t=%.1fs, ȷ=%d, ȷμ=%.1f, tReopt=%.1es, ReOptΣt=%.1fs\n"::Cstring;
        te::Cdouble,
        ȷ::Cint,
        ȷμ::Cdouble,
        tReopt::Cdouble,
        ReOptΣt::Cdouble
    )::Cint)
    @. Xl2 = round.(Xl; digits=2)
    println(Xl2)
end
_81(te) = rand() < 0.1 && @ccall(printf("t=%.1fs\n"::Cstring; te::Cdouble)::Cint)
function train_asyncly_shortly(mst, sub, WC_TIME)
    (; o, S, Xl, N, Θ, Vv, Σt), ch = mst, Channel{Int}(); CanSpawn = trues(S)
    Threads.@threads for n=sub
        Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(n.o), "MIPGap", TrainMIPGap)    
        Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(n.o), "TimeLimit", TrainTimeLim)
    end
    Settings.opt_ass_opt(mst, "prepare for async train"); load_che(mst)
    ȷ = nTook = 0; i = Wcpu; VΣ = sum(Vv); Vμ = VΣ/S
    for s=1:Wcpu _32(ch, s, sub, Xl, CanSpawn, Θ[s], Vμ) end
    t0 = time_ns()
    while true
        te = 1e-9(time_ns()-t0)
        if nTook > .5S
            ȷ > 0 ? _8(mst, Σt, ȷ, te) : _81(te)
            ȷ = nTook = 0
        elseif isready(ch) === false && ȷ > 0
            _8(mst, Σt, ȷ, te)
            ȷ = 0
        end
        te > WC_TIME && break
        s = take!(ch); nTook += 1; (;Ci,Cd) = n = sub[s]; CanSpawn[s] = true
        vio = Cd[N+1]
        if vio>Vμ
            gap = Settings.getmodeldblattr(n, "MIPGap")
            @ccall(printf("i=%lld | took s=%d, gap=%.1e, vio=%.4e, Vμ=%.4e\n"::Cstring;
                i::Clonglong,
                s::Cint,
                gap::Cdouble,
                vio::Cdouble, 
                Vμ::Cdouble
            )::Cint)
            Cd[N+1]=-1;Gurobi.GRBaddconstr(o,N+1,Ci,Cd,Cchar('<'),-Cd[end],C_NULL);ȷ+=1
        end
        u = vio>ϵ ? .3vio+.7ϵ : ϵ; VΣ += u-first(Vv); Vμ = VΣ/S; push!(Vv, u)
        while true
            i += 1; s = Settings.to_1S(i, S)
            CanSpawn[s] && (_32(ch, s, sub, Xl, CanSpawn, Θ[s], Vμ); break)
        end
    end
    for _=1:Wcpu
        s = take!(ch)
        (; Ci, Cd) = sub[s]
        vio = Cd[N+1]
        vio>Vμ && (Cd[N+1]=-1;Gurobi.GRBaddconstr(o,N+1,Ci,Cd,Cchar('<'),-Cd[end],C_NULL);ȷ+=1)
        u = vio>ϵ ? .5vio+.5ϵ : ϵ; VΣ += u-first(Vv); Vμ = VΣ/S; push!(Vv, u)
    end
end

function _32(ch, s, sub, mstXl, CanSpawn, Θs_che, Ev)
    sub[s].x_che .= mstXl
    CanSpawn[s], _ = false, Threads.@spawn(_31(ch, s, sub, Θs_che, Ev))
end
_31(ch, s, sub, Θs_che, Ev) = (_3(sub[s], Θs_che, Ev); put!(ch, s))
function _3(n, Θs_che, Ev) # ✅ SB_core
    (; o, N, Cd, Xl, Xl2, Bv, x_che) = n # local x_che should be loaded in advance
    Gurobi.GRBsetdblattrarray(o, "LB", 0, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 0, N, x_che)
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv), Bv)
    Settings.opt_ass_opt(n, "LP_relaxed_subproblem")
    Obj_LP = Settings.getmodeldblattr(n, "ObjVal")
    Pi_ast = Xl; Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Pi_ast)
    innprod = Pi_ast'x_che
    vio_LP = Obj_LP - Θs_che
    if vio_LP > Ev
        Cd[N+1] = vio_LP
        Cd[end] = Cn_LP = Obj_LP - innprod # This value aims to fit the Gurobi.addConstr API
        # when the training goes deeper, it's rare to still see separating LP cuts, so the following is rare
        @ccall(printf("¶ vio_LP = %.4e ¶\n"::Cstring; vio_LP::Cdouble)::Cint)
    else
        Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl2)
        Xl2 .= 1; Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl2)
        @. Xl2 = 0. - Pi_ast; Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
        Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv), Bv)
        Settings.opt_ass_time(n, "biased_MIP_subproblem")
        Cd[end] = Cn_MIP = Settings.getmodeldblattr(n, "ObjBound")
        Obj_MIP = Cn_MIP + innprod # This value facilitates comparison to the LP counterpart
        Cd[N+1] = vio_MIP = Obj_MIP - Θs_che
        # @ccall(printf("¶ Θs=%.4e, ObjLP=%.4e, ObjMIP=%.4e ¶\n"::Cstring;
        #     Θs_che::Cdouble, Obj_LP::Cdouble, Obj_MIP::Cdouble
        # )::Cint)
    end
    Cd[1:N] .= Pi_ast
end

function _4(ch, s, sub, x_che, VType) # ✅ eval_core
    (; o, N, Cd, Xl, Xl2, Bv) = n = sub[s]
    Gurobi.GRBsetdblattrarray(o, "LB", 0, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 0, N, x_che)
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
    Bv .= VType; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv), Bv)
    if Cchar(VType) === Cchar('C')
        Settings.opt_ass_opt(n, "eval LP subproblem")
        Cd[end] = ub_s = Settings.getmodeldblattr(n, "ObjVal")
        Cd[N+1] = 0.
    else
        Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(o), "MIPGap", EvalMIPGap)    
        Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(o), "TimeLimit", EvalTimeLim)    
        Settings.opt_ass_time(n, "eval MIP subproblem")
        Cd[end] = ub_s = Settings.getmodeldblattr(n, "ObjVal")
        Cd[N+1] = ub_s * Settings.getmodeldblattr(n, "MIPGap")
        Gurobi.GRBreset(o, 0)
    end
    put!(ch, s)
end
function eval_at_sync_point(mst, sub, VType)
    (; Xl, S, N), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst, "prepare it for eval")
    load_che(mst)
    lb, ub, GRBGapAbs = Settings.getmodeldblattr(mst, "ObjVal"), 0., 0.
    for s=1:WcpuEval Threads.@spawn(_4(ch, s, sub, Xl, VType)) end
    k = WcpuEval
    while true
        k === S && break
        Cd = sub[take!(ch)].Cd
        ub += Cd[end]
        GRBGapAbs += Cd[N+1]
        s = k += 1
        Threads.@spawn(_4(ch, s, sub, Xl, VType))
    end
    for _=1:WcpuEval
        Cd = sub[take!(ch)].Cd
        ub += Cd[end]
        GRBGapAbs += Cd[N+1]
    end
    ub /= S; GRBGapAbs /= S; agap = ub - lb; rgap = agap / ub
    printstyled("VType=$VType, lb=$lb, agap=$agap (GRBMIPErr=$GRBGapAbs), rgap=$rgap, ub=$ub\n"; color = :magenta)
end


function load_che(mst)
    (; o, N, S, Xl, Θ) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 0, S,  Θ)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
end
set_mst_θ_obj(mst) = ((;o,S,Θ)=mst; Θ .= 1/S; Gurobi.GRBsetdblattrarray(o, "Obj", 0, S, Θ))
function _2(o, sub, ch, N)
    (; Ci, Cd) = sub[ take!(ch) ]
    Cd[N+1] = -1
    Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL)
end
function MIP_warm(mst, sub) # only use Classic Benders cuts to warm up
    (; o, S, Xl, N), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst, "very first time, mst")
    load_che(mst); Threads.@threads for s=eachindex(sub)
        sub[s].x_che .= Xl
    end # Since this is performance critical, allocate a Channel
    for s=1:Wcpu Threads.@spawn(_31(ch, s, sub, -Inf, -Inf)) end
    i = Wcpu
    while true
        i+1 > S && break
        _2(o, sub, ch, N)
        s = i+=1
        Threads.@spawn(_31(ch, s, sub, -Inf, -Inf))
    end
    for _=1:Wcpu _2(o, sub, ch, N) end
    set_mst_θ_obj(mst)
    Settings.opt_ass_opt(mst, "The 1st time optimize with Objective properly set")
    load_che(mst)
end

end
