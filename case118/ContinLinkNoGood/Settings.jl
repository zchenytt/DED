"""
21: deep blue
27: blue
36: green
208: brown
"""
module Settings
import JuMP, Gurobi

# ✅ column generation
# Gurobi.GRBchgcoeffs(o, 2, Ci_con_vec, Ci_var_vec, Cd)
# Settings.setxdblattrelement(mst, N, "RHS", 0.)
# Gurobi.GRBsetdblattrarray(o, "RHS", start, len, Cd)
# Gurobi.GRBaddvar(o, 1+1, Cint[0,1], Cd, Qy_col, 0., 1e100, Cchar('C'), C_NULL)

# ✅ Suboptimal
# Settings.getmodelintattr(m, "SolCount")
# Gurobi.GRBsetintparam(Gurobi.GRBgetenv(m.o), "SolutionNumber", 0)
# "X" -> "PoolNX"

# ✅ param
# Gurobi.GRBsetintparam(Gurobi.GRBgetenv(mst.o), "OutputFlag", 1)
# Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(mst.o), "TimeLimit", 20.4)
# Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(mst.o), "NodeLimit", 0 or 1e100) # 0 is compute only at root node
# Gurobi.GRBsetdblparam(Gurobi.GRBgetenv(mst.o), "MIPGap", 0) # This can be written as 0::Int

"""
create Envs
"""
const C = Dict{String,Any}("Threads"=>1,"OutputFlag"=>0)
Env() = Gurobi.Env(C)
function Env(N::Int)
    v = Vector{Gurobi.Env}(undef, N)
    Threads.@threads for i=eachindex(v)
        v[i] = Env()
    end
    v
end

"""
create Models
"""
function Model(e::Gurobi.Env)
    m = JuMP.direct_model(Gurobi.Optimizer(e))
    JuMP.set_string_names_on_creation(m, false)
    m
end
function Model!(mv::Vector{JuMP.Model}, i, #=with existing ones=# ev::Vector{Gurobi.Env})
    m = JuMP.direct_model(Gurobi.Optimizer(ev[i]))
    JuMP.set_string_names_on_creation(m, false)
    mv[i] = m
end
Model!(mv::Vector{JuMP.Model}, #=with existing ones=# ev::Vector{Gurobi.Env}) = Threads.@threads for i=eachindex(mv)
    Model!(mv, i, ev)
end

# ✅ e.g. setting `x::JuMP.VariableRef`'s "Obj" Attribute
getxdblattrelement(m, i::Integer, str) = (r=m.refd;Gurobi.GRBgetdblattrelement(m.o, str,i,r);r.x)
setxdblattrelement(m, i::Integer, str, v) = Gurobi.GRBsetdblattrelement(m.o, str, i, v)
setxcharattrelement(m, i::Integer, str, v) = Gurobi.GRBsetcharattrelement(m.o, str, i, v) # e.g. Cchar('B')
# GRBgetcharattrelement(m.o, "VType", 0, &first_one);

# ✅ e.g. add a row
# Gurobi.GRBaddconstr(o, read_len, ci, cd, Cchar('=' or '>'), Gn, C_NULL)

# ✅ e.g. query "X" of the trial vec in 1st-stage, and then fix the copy vec in 2nd-stage
# Gurobi.GRBgetdblattrarray(o, "X", start, len, Ptr)
# Gurobi.GRBsetdblattrarray(o, "LB", start, len, Ptr)
# Gurobi.GRBsetcharattrarray(o, "VType", start, len, fill(Cchar('B'), len)) # ✅ e.g. set variables to binary

# ✅ e.g. query ObjVal/Runtime
function getmodeldblattr(m, str)
    r = m.refd
    Gurobi.GRBgetdblattr(m.o, str, r)
    r.x
end
# ✅ e.g. query NumConstrs/NumVars
function getmodelintattr(m, str)
    r = m.refi
    Gurobi.GRBgetintattr(m.o, str, r)
    r.x
end
opt_and_ter(m) = ((; o, refi)=m;Gurobi.GRBoptimize(o);Gurobi.GRBgetintattr(o, "Status", refi); refi.x)
opt_ass_time(m) = opt_ass_time(m, "")
opt_ass_opt(m) = opt_ass_opt(m, "")
function opt_ass_opt(m, str)
    t=opt_and_ter(m)
    t == 2 || (@error(string(str, "> ter = ", t)); error(t))
end
function opt_ass_time(m, str)
    t=opt_and_ter(m)
    (t == 2 || t == 9) || (@error(string(str, "> ter = ", t)); error(t))
end

to_1S(i, S) = (i-1) % S + 1 # this projects any positive integer to 1:S
isfrac(x) = ≈(x, 0.5; atol = 0.49999)
function printinfo()
    a,b,c,d = Threads.nthreads(:default),Threads.nthreads(:interactive),Threads.threadid(),Threads.threadpool()
    println("Settings> Threads=($a,$b), tid=$c, pool=$d")
end
_gcc(o,x) = Gurobi.c_column(o,JuMP.index(x))
end

# ✅ Status: 3==INFEASIBLE, 5==UNBOUNDED, 4==(3||5), 9==TimeLimit, 1==Not_solved_yet
