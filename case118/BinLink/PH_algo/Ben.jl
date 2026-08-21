module Ben
import ..Models, ..Settings, Gurobi

const EVAL_MIPGap = 1e-4
const EVAL_TimLim = 45.
const MAIN_MIPGap = 2e-4
const MAIN_TimLim = 15. # 🟥 if MIP-sub is too hard, tune its hardness down!
const W_fg, W_bg = Threads.nthreads(:interactive), Threads.nthreads(:default)

function para_build(ch, envs, a...)
    i = 1; S = length(envs)
    @assert W_bg ≤ S
    for _=1:W_bg Models._1(ch, i, envs, a...); i+=1 end
    while i ≤ S
        take!(ch)
        Models._1(ch, i, envs, a...)
        i += 1
    end
    for _=1:W_bg take!(ch) end
end

function _6(ch, sub, s, Xl)
    (; o, N, Bv, Pi, Xl2, Cd) = n = sub[s]
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
    Xl2 .= 0;  Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
               Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
               Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
    Settings.opt_ass_opt(n, "no bias sub-LP")
    Cd[N+1] = obj = Settings.getmodeldblattr(n, "ObjVal")
    Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Pi)
    Cd[1:N] .= Pi
    Cd[end] = Cn = obj - Pi'Xl
    put!(ch, s)
end
_7(ch, sub, s, Xl) = Threads.@spawn(_6(ch, sub, s, Xl))
function add_a_single_θ_cut(mst, sub; mstVType)
    (; o, S, N, Θ, Xl, Σt, Bv, ch, Ci, Cd) = mst
    mstVType === 'B' || error("not planned")
    Bv .= mstVType; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv)
    Settings.opt_ass_opt(mst)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
    @. Xl = round(Xl); xHa = hash(Xl)
    θμ_che = Settings.getxdblattrelement(mst, S+N, "X")
    i = 1
    for _=1:W_bg _7(ch, sub, i, Xl); i+=1 end
    while i ≤ S
        take!(ch)
        _7(ch, sub, i, Xl); i+=1
    end
    for _=1:W_bg take!(ch) end
    ub_LP = sum(n.Cd[N+1] for n=sub)/S
    ν = ub_LP - θμ_che
    @ccall(printf("x_int_che = %llx, vio_by_LP = %.4e\n"::Cstring; xHa::UInt64, ν::Cdouble)::Cint)
    cut_added = ν > ϵ
    if cut_added
        Cd .= 0
        for n=sub Cd .+= n.Cd end
        Cd ./= S
        Cd[N+1] = -1
        Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL)
    end
    cut_added, xHa
end

function get_ub(s::Int, Xl::Vector{Float64}, sub) # ub of Q_s(x)
    (; o, N, Bv, Pi, Xl2, Cd) = n = sub[s]
    Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(o), "MIPGap", EVAL_MIPGap)
    Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(o), "TimeLimit", EVAL_TimLim)
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 0, N, Bv)
    Xl2 .= 0;  Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
               Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
               Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv)-N, Bv)
    Settings.opt_ass_time(n, "no bias sub-MIP eval")
    Cd[end] = ub = Settings.getmodeldblattr(n, "ObjVal")
    Cd[N+1] = agap = Settings.getmodeldblattr(n, "MIPGap")abs(ub)
end
_9(s, Xl, sub, ch) = (get_ub(s, Xl, sub); Gurobi.GRBreset(sub[s].o, 0); put!(ch, s))
_8(s, Xl, sub, ch) = Threads.@spawn(_9(s, Xl, sub, ch))

function get_ub(Xl::Vector{Float64}, sub, ch)
    i, S = 1, length(sub)
    for _=1:W_bg _8(i, Xl, sub, ch); i+=1 end
    while i ≤ S
        take!(ch)
        _8(i, Xl, sub, ch)
        i+=1
    end
    for _=1:W_bg take!(ch) end
    ub = sum(n.Cd[end] for n=sub)/length(sub)
    errsum = sum(n.Cd[end-1] for n=sub)
    erravg = errsum/length(sub)
    printstyled("ub = $ub, erravg = $erravg, errsum = $errsum\n"; color = 27)
    ub
end

function _4(ch, s, Θs, n, Vμ)
    (; o, N, Bv, Pi, Xl, Xl2, Cd, rtCnt) = n
    _en2 = Gurobi.GRBgetenv(o)
    Gurobi.GRBsetdblparam(_en2, "TimeLimit", MAIN_TimLim)
    Gurobi.GRBsetdblparam(_en2, "MIPGap", MAIN_MIPGap)
    Bv .= 'C';       Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
    Xl2 .= 0;        Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
                    Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
                    Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
    Settings.opt_ass_opt(n, "no bias - subLP")
    Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Pi)
    Cd[1:N] .= Pi
    @. Xl2 = 0 - Pi; Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
    Xl2 .= 0;        Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl2)
    Xl2 .= 1;        Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl2)
    Bv .= 'B';       Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
    Settings.opt_ass_time(n, "biased - subMIP")
    GRBrgap, Rt = Settings.getmodeldblattr(n, "MIPGap"), Settings.getmodeldblattr(n, "Runtime")
    Rt+0.1>MAIN_TimLim ? rtCnt.o += 1 : rtCnt.n += 1
    GRBagap = abs(Settings.getmodeldblattr(n, "ObjVal"))GRBrgap
    Cd[end] = Obn = Settings.getmodeldblattr(n, "ObjBound")
    Cd[N+1] = ν = Pi'Xl + Obn - Θs # level of separation, at the current x_che
    (rand()<5e-5 && ν>0) && @ccall(printf("  GRBgap=(r=%.1e, a=%.1f; t=%.1f), vio=%.4e <s=%d\n"::Cstring;GRBrgap::Cdouble,GRBagap::Cdouble,Rt::Cdouble,ν::Cdouble,s::Cint)::Cint)
    put!(ch, s)
end
function _5(ch, s, sub, mst, CanSpawn, Vμ)
    Θs, n = mst.Θ[s], sub[s]
    n.Xl .= mst.Xl
    CanSpawn[s], _ = false, Threads.@spawn(_4(ch, s, Θs, n, Vμ))
end
function root_train(mst, sub, SeqSeconds, AbsSeconds)
    (; o, S, N, Θ, Σt, Bv, Vv, CanSpawn, ch, rng, proceed) = mst
    CanSpawn .= true;       VΣ=sum(Vv); Vμ=VΣ/S; α = .5
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv)
    Settings.opt_ass_opt(mst); load_che(mst); lb = Settings.getmodeldblattr(mst, "ObjBound")
    ȷ, t0 = 0, time_ns(); i = i0 = rand(rng, 1:S)
    for _=1:W_bg # load the hardware
        _5(ch, i, sub, mst, CanSpawn, Vμ)
        i = ifelse(i<S, i+1, 1)
    end
    iA = iB = ter = 0
    while true # proceed.value
        (t = 1e-9(time_ns() - t0)) > AbsSeconds && (ter = 3; break)
        if ȷ > 0 && (!isready(ch) || i === i0)
            Settings.opt_ass_opt(mst); load_che(mst)
            lb = Settings.getmodeldblattr(mst, "ObjBound")
            (st = Σt.x += Settings.getmodeldblattr(mst, "Runtime")) > SeqSeconds && (ter = 4; break)
            nCμ = Settings.getmodelintattr(mst, "NumConstrs")/S
            rand()<0.0001 && @ccall(printf("M> t=%.1fs, ȷ=%d, nCμ=%.1f, st=%.1fs, lb=%.1f, Vμ=%.1e\n"::Cstring;t::Cdouble,ȷ::Cint,nCμ::Cdouble,st::Cdouble,lb::Cdouble,Vμ::Cdouble)::Cint)
            ȷ = iA = 0
        else
            if iA === S
                @ccall(printf("idle-M> t=%.1fs, lb=%.1f, Vμ=%.1e\n"::Cstring;t::Cdouble,lb::Cdouble,Vμ::Cdouble)::Cint)
                iB += 1
                iB === 2 && (ter = 2; break)
                iA = 0
            else
                iA += 1
            end
        end
        s = take!(ch)
        (; Ci, Cd), CanSpawn[s] = sub[s], true
        ν, Cd[N+1] = Cd[N+1], -1.
        ν>Vμ && (Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL); ȷ+=1; iB=0)
        ν = ν>ϵ ? (α)ν+(1-α)ϵ : ϵ; VΣ=VΣ-first(Vv)+ν; Vμ=VΣ/S; push!(Vv,ν)
        while true
            brkFlag = CanSpawn[i]
            brkFlag && _5(ch, i, sub, mst, CanSpawn, Vμ)
            i = ifelse(i<S, i+1, 1)
            brkFlag && break
        end
    end
    for _=1:W_bg # wait for all threads to finish
        s = take!(ch)
        (; Ci, Cd), CanSpawn[s] = sub[s], true
        ν, Cd[N+1] = Cd[N+1], -1.
        ν>Vμ && (Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL); ȷ+=1; iB=0)
        ν = ν>ϵ ? (α)ν+(1-α)ϵ : ϵ; VΣ=VΣ-first(Vv)+ν; Vμ=VΣ/S; push!(Vv,ν)
    end
    nRt, oRt = sum(n.rtCnt.n for n=sub), sum(n.rtCnt.o for n=sub)
    ħ = round(oRt/(oRt+nRt); digits = 2)
    printstyled("root_train> ter=$ter, seqTime=$(Σt.x), Hardness=$ħ\n";color=208)
end
function run_a_round(mst, sub; mstVType, addcut)
    (; o, S, N, Θ, Bv, CanSpawn, ch, Xl) = mst
    if mstVType === 'B'
        Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "OutputFlag", 1)
    elseif mstVType === 'C'
        Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "OutputFlag", 0)
    else
        error()
    end
    Bv .= mstVType; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv)
    Settings.opt_ass_opt(mst); load_che(mst)
    if addcut
        i = 1
        for _=1:W_bg
            _5(ch, i, sub, mst, CanSpawn, 0.)
            i += 1
        end
        while i ≤ S
            (; Ci, Cd) = sub[take!(ch)]
            Cd[N+1] = -1; Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL)
            _5(ch, i, sub, mst, CanSpawn, 0.)
            i += 1
        end
        for _=1:W_bg
            (; Ci, Cd) = sub[take!(ch)]
            Cd[N+1] = -1; Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL)
        end
        Θ .= 1/S; Gurobi.GRBsetdblattrarray(o, "Obj", 0, S, Θ)
        Settings.opt_ass_opt(mst); load_che(mst)
    end
    mstVType === 'B' && @. Xl = round(Xl)
    lb = Settings.getmodeldblattr(mst, "ObjBound")
    printstyled("mst's lb = $lb, xHa = $(repr(hash(Xl)))\n"; color = 30)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "OutputFlag", 0)
    lb
end

function leaf_train(mst, sub)
    Models.retrofit!(mst)
    h_old = UInt64(0)
    for k = 1:typemax(Int)
        cut_added, h = add_a_single_θ_cut(mst, sub; mstVType='B')
        cut_added && continue
        h_old === h && break
        h_old = h
    end
end


function load_che(mst)
    (; o, S, N, Θ, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 0, S,  Θ)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
end

end
