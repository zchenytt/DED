"""
Approximate Footprint: It spans a region of roughly 500 km in length by 350 km in width (approximately 175,000 km^2).

original generators capacity_sum = 65.15
original load = 42.42
added wind capacity = 40.2 (wind penetration = 38%)
Thereby we can scale up the load to 1.94x (We've verified that the transmission capacity is sufficient after scaling up)

Due to these modifications, RateA needs to be enlarged properly.

C_energy:C_reserve_1st:C_windcur:C_shed = 1:5:20:100  (large C_reserve_1st is due to opportunity cost of withholding generation capacity)
"""
module Case118
import PowerModels, PGLib
const CaseName = "pglib_opf_case118_ieee.m"
const B, N, L = 186, 118, 99
const Load_Multiplier = 1.7 # to scale up the load to match the new wind capacity
const RateA_Multiplier = 1.0
const LookAheadHours = 4 # We always look ahead 4 hours; Then we have 2 possibilities: (5min)48 _or_ (15min)16

# Nominal load; Time window: 16:00 → 20:00 (17 anchors, 15-min interval)
const l16 = (
    [0.92, 0.92, 0.91, 0.91, 0.9, 0.9, 0.89, 0.89, 0.88, 0.88, 0.87, 0.87, 0.86, 0.86, 0.85, 0.85, 0.84],
    [1.0, 0.98, 0.96, 0.93, 0.9, 0.86, 0.82, 0.77, 0.72, 0.67, 0.62, 0.58, 0.54, 0.5, 0.47, 0.44, 0.42],
    [0.55, 0.6, 0.66, 0.73, 0.8, 0.87, 0.93, 0.97, 1.0, 0.99, 0.97, 0.94, 0.9, 0.85, 0.8, 0.75, 0.7]
);
get_load_type() = rand(eachindex(l16), L)
get_load(CaD) = [Load_Multiplier * CaD["load"][string(i)]["pd"] for i=1:L]
get_loadnode(CaD) = [CaD["load"][string(i)]["load_bus"] for i=1:L]

get_Case_Dict()=(PowerModels.silence();PowerModels.make_basic_network(PGLib.pglib(CaseName)))
_0!(F; a=9e-7)=for (i,e)=enumerate(F) if abs(e)<a setindex!(F, 0., i) end end
get_PTDF(CaD)=(F=PowerModels.calc_basic_ptdf_matrix(CaD);_0!(F);F)
get_ratea(CaD)= [RateA_Multiplier * CaD["branch"][string(i)]["rate_a"] for i = 1:B]

"""
New wind units are added to bus 8, 65, 116.
Their wind pMax has the ratio 23.90:24.51:72.18 (from the transmission line capacity at that node)
"""
get_windnode() = [8, 65, 116]
get_windPmax() = [23.9, 24.51, 72.18] / 3

"""
New storage units are added to bus 59 (load side) and bus 116 (source side)
"""
get_esnode() = [59,  116]
get_esEmax() = [34., 36.] # Assume they (fastest possible) take 2 hours to be fully charged
get_esPmax() = [17., 18.] # based on 1 hour
get_esPmin() = [5.,   5.] # based on 1 hour

"""
Only pmax>0 are dispatchable generators, we have 19 such generators
However, there's a UC upper layer so that not all of them are online-dispatchable during intraday time
to this end, we need a Agi::Vector{Int}, so that `[GPmax[g] for (i,g)=enumerate(Agi)]`
"""
get_Ggvec(CaD) = [g for g = 1:54 if CaD["gen"][string(g)]["pmax"] > 0]
get_Gnode(CaD, Ggvec) = [CaD["gen"][string(g)]["gen_bus"] for g = Ggvec]
get_GPmax(CaD, Ggvec) = [CaD["gen"][string(g)]["pmax"] for g = Ggvec]
get_GClin(CaD, Ggvec) = [0.001CaD["gen"][string(g)]["cost"][2] for g = Ggvec]

"This is a helper function"
function calculate_node_max_trans_cap(CaD, n) # e.g. to verify that each load can be served via transmission
    c = 0.
    for i = 1:B
        b = CaD["branch"][string(i)]
        if b["f_bus"] == n || b["t_bus"] == n
            c += b["rate_a"]
        end
    end
    c
end

end
