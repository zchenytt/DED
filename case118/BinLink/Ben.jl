module Ben
import ..Models, ..Settings, Gurobi

const W_naive = Sys.CPU_THREADS
const W_fg = 1
const W_bg = W_naive - W_fg
const ϵ = Models.SB_CUT_COT

# 🟥 The final global optimality gap relies on the hardness of the subMIPs that are solved to suboptimality
# if the subMIPs are coarsely solved, the final global optimality gap can be large (e.g. rgap = 0.24).
config_sub_mips(sub; MIPGap=9e-4, TimeLimit=12.) = Threads.@threads(for s=eachindex(sub)
    _en2 = Gurobi.GRBgetenv(sub[s].o)
    Gurobi.GRBsetdblparam(_en2, "MIPGap", MIPGap)
    Gurobi.GRBsetdblparam(_en2, "TimeLimit", TimeLimit)
end)

function MIP_warm(mst, sub)
    (; S, N, Θ, Σt) = mst
    @assert S > W_bg
    Settings.opt_ass_opt(mst); load_che(mst)
    config_sub_mips(sub)
    Threads.@threads for s = eachindex(sub)
        n = sub[s]
        n.Xl .= mst.Xl # make a local static copy
        (; o, Bv, Pi, Xl, Xl2, Ci, Cd) = n
        Bv .= 'C';          Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
        Xl2 .= 0;           Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
                            Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
                            Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
        Settings.opt_ass_opt(n, "warm - LP")
        Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Pi)
        @. Xl2 = 0 - Pi;    Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
        Xl2 .= 0;           Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl2)
        Xl2 .= 1;           Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl2)
        Bv .= 'B';          Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
        Settings.opt_ass_time(n, "warm - biased subMIP")
        Cd[1:N] .= Pi
        Cd[end]  = Obn = Settings.getmodeldblattr(n, "ObjBound")
        Cd[N+1] = -1
    end
    for n=sub Gurobi.GRBaddconstr(mst.o, N+1, n.Ci, n.Cd, Cchar('<'), -n.Cd[end], C_NULL) end
    Θ .= 1/S; Gurobi.GRBsetdblattrarray(mst.o, "Obj", 0, S, Θ)
    Settings.opt_ass_opt(mst, "The 1st time optimize with Objective properly set")
    load_che(mst)
    Σt.x += Settings.getmodeldblattr(mst, "Runtime")
end
function get_lb_ub(mst, sub)
    (; S, N, Θ) = mst
    mst.Bv .= 'B'; Gurobi.GRBsetcharattrarray(mst.o, "VType", S, N, mst.Bv)
    Settings.opt_ass_opt(mst); load_che(mst)
    config_sub_mips(sub; MIPGap=1e-4, TimeLimit=50.)
    Threads.@threads for s = eachindex(sub)
        (; o, Bv, Pi, Xl, Xl2, Ci, Cd) = n = sub[s]
        @. Xl = round(mst.Xl)
        Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 0, N, Bv)
        Xl2 .= 0;  Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
                   Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
                   Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
        Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", N, length(Bv)-N, Bv)
        Settings.opt_ass_time(n, "no bias subMIP")
    end
    GRBagapErr = ub = 0.
    for (s,n)=enumerate(sub)
        ub += ObjVal = Settings.getmodeldblattr(n, "ObjVal")
        GRBagapErr += Settings.getmodeldblattr(n, "MIPGap")abs(ObjVal)
        Θ[s] = ObjVal - Θ[s] # update the agap
    end
    (agapmin, agapmax), agapsum = extrema(Θ), sum(Θ)
    rgap, agap = agapsum/ub, agapsum/S
    ub /= S
    printstyled(
        "ub = $ub, rgap = $rgap; agap=(($agapmin,$agapmax), ($agapsum - $GRBagapErr)/$S)",
        " for other details, check `mst.Θ`.\n"; color = :magenta
    )
end

function _4(ch, s, Θs, n, Vμ)
    (; o, N, Bv, Pi, Xl, Xl2, Cd) = n
    Bv .= 'C';          Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
    Xl2 .= 0;           Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
    Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl)
    Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl)
    Settings.opt_ass_opt(n, "train - subLP")
    Gurobi.GRBgetdblattrarray(o, "RC", 0, N, Pi)
    @. Xl2 = 0 - Pi;    Gurobi.GRBsetdblattrarray(o, "Obj", 0, N, Xl2)
    Xl2 .= 0;           Gurobi.GRBsetdblattrarray(o, "LB", 0, N, Xl2)
    Xl2 .= 1;           Gurobi.GRBsetdblattrarray(o, "UB", 0, N, Xl2)
    Bv .= 'B';          Gurobi.GRBsetcharattrarray(o, "VType", 0, length(Bv), Bv)
    Settings.opt_ass_time(n, "train - subMIP")
    Cd[1:N] .= Pi
    Cd[end]  = Obn = Settings.getmodeldblattr(n, "ObjBound")
    GRBrgap = Settings.getmodeldblattr(n, "MIPGap")
    GRBagap = abs(Settings.getmodeldblattr(n, "ObjVal"))GRBrgap
    Rt = Settings.getmodeldblattr(n, "Runtime")
    Cd[N+1]  = vio = Pi'Xl + Obn - Θs
    vio > Vμ/2 && @ccall(printf("s=%d, GRBgap=(r=%.1e, a=%.1f; t=%.1f), vio=%.4e ◀\n"::Cstring;
        s::Cint,
        GRBrgap::Cdouble,
        GRBagap::Cdouble,
        Rt::Cdouble,
        vio::Cdouble
    )::Cint)
    put!(ch, s)
end
function _5(ch, s, sub, mst, CanSpawn, Vμ)
    Θs, n = mst.Θ[s], sub[s]
    n.Xl .= mst.Xl
    CanSpawn[s], _ = false, Threads.@spawn(_4(ch, s, Θs, n, Vμ))
end
function main(mst, sub, SeqSeconds, AbsSeconds)
    (; o, S, N, Θ, Σt, Bv, Vv) = mst
    VΣ=sum(Vv); Vμ=VΣ/S; α = .5
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv)
    config_sub_mips(sub)
    ch, CanSpawn = Channel{Int}(), trues(S)
    Settings.opt_ass_opt(mst); load_che(mst)  
    i, spnˈd, ȷ = rand(1:S), 0, 0
    while true # load the hardware
        spnˈd == W_bg && break
        _5(ch, i, sub, mst, CanSpawn, Vμ); spnˈd += 1
        i = ifelse(i<S, i+1, 1)
    end
    spnˈold = spnˈd; t0 = time_ns()
    while true
        (t = 1e-9(time_ns() - t0)) > AbsSeconds && break
        if ȷ > 0 && (isready(ch) === false || spnˈd > spnˈold + W_bg)
            Settings.opt_ass_opt(mst, "mst_reopt"); load_che(mst)
            st = Σt.x += Settings.getmodeldblattr(mst, "Runtime")
            st > SeqSeconds && break
            nCμ = Settings.getmodelintattr(mst, "NumConstrs")/S
            @ccall(printf("▶ t=%.1fs, ȷ=%d, nCμ=%.1f, st=%.1fs, Vμ=%.1e\n"::Cstring;
            t::Cdouble, ȷ::Cint, nCμ::Cdouble, st::Cdouble, Vμ::Cdouble)::Cint)
            ȷ = 0
            spnˈold = spnˈd
        end
        s = take!(ch)
        CanSpawn[s]=true
        (; Ci, Cd) = sub[s]
        ν, Cd[N+1] = Cd[N+1], -1.
        ν>Vμ && (Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL); ȷ += 1)
        ν = ν>ϵ ? (α)ν+(1-α)ϵ : ϵ; VΣ=VΣ-first(Vv)+ν; Vμ=VΣ/S; push!(Vv,ν)
        while true
            brkFlag = CanSpawn[i]
            brkFlag && (_5(ch, i, sub, mst, CanSpawn, Vμ); spnˈd += 1)
            i = ifelse(i<S, i+1, 1)
            brkFlag && break
        end
    end
    for _=1:W_bg # wait for all threads to finish
        s = take!(ch)
        (; Ci, Cd) = sub[s]
        ν, Cd[N+1] = Cd[N+1], -1.
        ν>Vμ && Gurobi.GRBaddconstr(o, N+1, Ci, Cd, Cchar('<'), -Cd[end], C_NULL)
        ν = ν>ϵ ? (α)ν+(1-α)ϵ : ϵ; VΣ=VΣ-first(Vv)+ν; Vμ=VΣ/S; push!(Vv,ν)
    end
end
function load_che(mst)
    (; o, S, N, Θ, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 0, S,  Θ)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
end

end
