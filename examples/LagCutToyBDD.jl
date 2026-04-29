include("src/Settings.jl");
import Gurobi, JuMP, Random; Random.seed!(2);
function build_master(m)
    JuMP.@variable(m, t)
    JuMP.@variable(m, 0 <= y <= 1)
    JuMP.@objective(m, Min, 0)
end;
function build_subproblem(m)
    JuMP.@variable(m, y, Bin)
    JuMP.@variable(m, x >= 0)
    JuMP.@constraint(m, x + 15y >= 8)
    JuMP.@constraint(m, 3x + 10y >= 13)
    JuMP.@constraint(m, x + 10y >= 7)
    JuMP.@constraint(m, 2x - 10y >= -1)
    JuMP.@constraint(m, 2x - 70y >= -49)
    JuMP.@objective(m, Min, x)
end;
function get_mst_sub()
    genv = Settings.Env();mst = Settings.Model(genv);o = mst.moi_backend
    build_master(mst)
    mst = (m = mst, o = o, refi = Ref{Cint}(-99), refd = Ref{Cdouble}(NaN))
    genv = Settings.Env();sub = Settings.Model(genv);o = sub.moi_backend
    build_subproblem(sub)
    sub = (m = sub, o = o, refi = Ref{Cint}(-99), refd = Ref{Cdouble}(NaN))
    genv = Settings.Env();pas = Settings.Model(genv);o = pas.moi_backend
    build_pas(pas)
    pas = (m = pas, o = o, refi = Ref{Cint}(-99), refd = Ref{Cdouble}(NaN))
    mst, sub, pas
end;
function build_pas(pas) # pi ascending; The hardest point is that it is parameterized by the master's trial point and it doesn't have an upper bound usually
    JuMP.@variable(pas, h) # hypograph variable
    JuMP.@variable(pas, p) # pi
    # JuMP.@constraint(m, h <= -0.0 * p + 8)
    # JuMP.@constraint(m, h <= -1.0 * p + 10.5)
    JuMP.@objective(pas, Max, h)
end;
function get_check_and_hat_initial(mst, sub)
    y_check = Settings.getxdblattrelement(mst, 1, "X")
    Settings.setxdblattrelement(sub, 0, "Obj", 0.)
    Settings.setxdblattrelement(sub, 0, "LB", y_check)
    Settings.setxdblattrelement(sub, 0, "UB", y_check)
    Settings.setxcharattrelement(sub, 0, "VType", Cchar('C'))
    Settings.opt_ass_opt(sub)
    p_hat = Settings.getxdblattrelement(sub, 0, "RC")
    y_check, p_hat
end;
function get_subLB_and_col(sub, p_hat)
    Settings.setxdblattrelement(sub, 0, "Obj", -1.0 * p_hat);
    Settings.setxdblattrelement(sub, 0, "LB", 0.0);
    Settings.setxdblattrelement(sub, 0, "UB", 1.0);
    Settings.setxcharattrelement(sub, 0, "VType", Cchar('B'));
    Settings.opt_ass_opt(sub);
    subLB = Settings.getmodeldblattr(sub, "ObjBound")
    y_col = Settings.getxdblattrelement(sub, 0, "X")
    Qy_col = Settings.getxdblattrelement(sub, 1, "X")
    subLB, y_col, Qy_col
end;

mst, sub, pas = get_mst_sub();
Settings.opt_ass_opt(mst);
# iter 1
y_check, p_hat = get_check_and_hat_initial(mst, sub) # (0.0, -15.0)
subLB, y_col, Qy_col = get_subLB_and_col(sub, p_hat) # (8.0, 0.0, 8.0)
JuMP.@constraint(pas.m, pas.m[:h] <= -y_col * pas.m[:p] + Qy_col); # column generation => dual row generation
JuMP.@constraint(mst.m, mst.m[:t] >= p_hat * mst.m[:y] + subLB); # add optimality cut to master
Settings.setxdblattrelement(mst, 0, "Obj", 1.0); # activate theta in master
Settings.opt_ass_opt(mst);
# iter 2
y_check, p_hat = get_check_and_hat_initial(mst, sub) # (1.0, 35.0)
subLB, y_col, Qy_col = get_subLB_and_col(sub, p_hat) # (-24.5, 1.0, 10.5)
JuMP.@constraint(pas.m, pas.m[:h] <= -y_col * pas.m[:p] + Qy_col); # column generation => dual row generation
JuMP.@constraint(mst.m, mst.m[:t] >= p_hat * mst.m[:y] + subLB); # add optimality cut to master
Settings.opt_ass_opt(mst);
# iter 3
y_check, p_hat = get_check_and_hat_initial(mst, sub)
Settings.setxdblattrelement(pas, 1, "Obj", y_check); # the pi-ascending problem is parameterized by check value
Settings.opt_ass_opt(pas);
JuMP.value(pas.m[:p]) # 2.5
