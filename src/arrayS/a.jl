
# TODO add line flow limits to the zones
# TODO there is another way of implementing: allocate subproblem models number fixed to the hardware (this is memory efficient)
# We need to compare these two situations
# see which one has better overall performance

import Random; Random.seed!(hash(273));
const t, T, S = 1, 24, 10000
Ths = min(S, 253)
MaxSec = 45*60

using Ded
Settings.printinfo()
import JuMP, Gurobi
const Other, Reserve, Demand, Load = Ded.LoadGen.get_case2383(T); # Demand.Zone = [4, 5, 6]
const Wind = Ded.WindGen.get_case2383(t, T, S);
const tks = similar(Wind.S, Task);
const Lk = (ref = ReentrantLock(), mst = ReentrantLock());
const mst = Settings.Model();
const inn = Settings.Model(tks);

Build._mst(mst, t, T, S, Other, Reserve, Wind, Demand, Load) # this goes before building inn
Build.inner_model(tks, inn, t, T, Reserve, Wind, Demand, Load, mst)

Build.add_pibound_cut(tks, mst, Lk.mst, t, T, Other, Reserve, Wind, Demand, Load)
Settings.opt_and_ter(mst) === JuMP.OPTIMAL || error()
const ref = Ref(Build.get_nt(mst));
println("Initial LB = $(ref.x.lb)")
const Ncuts = Threads.Atomic{Int}(0);
const otr = Ref{Task}();
const evt = Base.Event(true);

otr.x = Threads.@spawn(identity); wait(otr.x)
stk = Train.ss(evt, 1, otr, tks, Lk, ref, mst, inn, Ncuts, Ths, MaxSec)
# Training asyncly...
Train.wait3(stk, otr, tks)
Train.post_check(tks, inn, mst)
















function train(N, inn, a...)
    for k = 1:N
        println("k = $k")
        for s = 1:S
            Train.scene_subprocedure(inn, s, a...)
        end
        Train.mst_procedure(a...)
    end
end
train(40, inn, mst, Ncuts, Lk, ref)






