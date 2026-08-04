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
# Nominal load; Time window: 16:00 → 20:00 (17 anchors, 15-min interval)
const l16 = (
    [0.92, 0.92, 0.91, 0.91, 0.9, 0.9, 0.89, 0.89, 0.88, 0.88, 0.87, 0.87, 0.86, 0.86, 0.85, 0.85, 0.84],
    [1.0, 0.98, 0.96, 0.93, 0.9, 0.86, 0.82, 0.77, 0.72, 0.67, 0.62, 0.58, 0.54, 0.5, 0.47, 0.44, 0.42],
    [0.55, 0.6, 0.66, 0.73, 0.8, 0.87, 0.93, 0.97, 1.0, 0.99, 0.97, 0.94, 0.9, 0.85, 0.8, 0.75, 0.7]
);

const EVnode = [42, 54, 59, 80, 90, 116]
# Physically according to the formula, EVLoad should be _energy_,
# In JuMP, for visual clarity, we equivalently use discharging _power_
const EVLmax = 3.405 * [0.96, 1.13, 2.77, 1.3, 1.63, 1.84]
const EVL01 = [l16..., map(reverse, l16)...]
const EVEmax = [2.85, 4.24, 2.82, 1.13, 3.31, 2.52]
const EVEini = [1.6815, 2.5016, 1.6638, 0.6667, 1.9529, 1.0416]
const Bxini = Bool[1, 1, 0, 1, 1, 0]

const EEini = [5.9, 14.16, 29.5, 9.84]

const B, N, L = 186, 118, 99-length(EVnode)
const Load_Multiplier = 1.0 # to scale up the load to match the new wind capacity
const RateA_Multiplier = 1.0
const LookAheadHours = 4 # We always look ahead 4 hours; Then we have 2 possibilities: (5min)48 _or_ (15min)16
const MyRateA = [0.045, 0.474, 1.4, 0.889, 1.277, 0.809, 0.007, 4.789, 0.007, 1.065, 1.234, 1.009, 0.183, 0.074, 0.643, 0.667, 0.534, 0.361, 0.408, 0.468, 1.031, 0.268, 0.897, 0.363, 0.29, 0.53, 0.117, 0.077, 0.169, 0.263, 0.489, 1.022, 0.597, 0.095, 0.122, 2.572, 4.614, 1.022, 0.718, 0.36, 0.445, 0.061, 0.202, 0.322, 0.198, 0.097, 0.334, 0.117, 0.27, 0.982, 3.383, 1.176, 1.047, 1.083, 0.935, 0.861, 0.716, 0.668, 0.162, 0.337, 0.048, 0.254, 0.408, 0.09, 0.538, 0.801, 0.801, 0.272, 0.266, 0.743, 0.922, 0.396, 0.244, 0.033, 0.682, 0.662, 0.109, 0.843, 0.219, 0.49, 0.6, 0.314, 0.387, 0.422, 0.336, 0.351, 0.426, 1.245, 1.337, 1.632, 0.377, 0.163, 4.881, 4.881, 3.002, 2.688, 7.86, 2.312, 2.312, 0.492, 0.374, 5.894, 0.626, 7.163, 1.166, 1.09, 7.752, 1.848, 0.231, 0.327, 0.151, 0.269, 0.06, 0.476, 0.396, 2.241, 0.251, 0.36, 4.249, 0.822, 1.031, 0.335, 0.423, 0.204, 0.399, 10.96, 10.96, 3.367, 2.4, 0.911, 1.314, 0.814, 0.194, 0.007, 0.957, 0.783, 0.585, 0.954, 1.831, 1.074, 1.154, 0.361, 1.122, 1.077, 1.11, 1.164, 1.458, 1.753, 0.51, 1.831, 1.827, 2.027, 1.942, 0.496, 0.861, 1.777, 1.685, 1.724, 1.57, 0.629, 0.374, 0.416, 1.582, 0.635, 0.262, 0.376, 0.647, 0.561, 0.008, 0.23, 0.404, 0.232, 0.386, 0.743, 0.328, 0.007, 0.68, 0.54, 0.496, 0.236, 0.098, 0.156, 25.244, 0.184, 1.233, 0.965]



get_load_type() = rand(eachindex(l16), L)
get_load(CaD) = [Load_Multiplier * CaD["load"][string(i)]["pd"] for i=1:99 if CaD["load"][string(i)]["load_bus"] ∉ EVnode]
get_loadnode(CaD) = [CaD["load"][string(i)]["load_bus"] for i=1:99 if CaD["load"][string(i)]["load_bus"] ∉ EVnode]

get_Case_Dict()=(PowerModels.silence();PowerModels.make_basic_network(PGLib.pglib(CaseName)))
_0!(F; a=9e-7)=for (i,e)=enumerate(F) if abs(e)<a setindex!(F, 0., i) end end
get_PTDF(CaD)=(F=PowerModels.calc_basic_ptdf_matrix(CaD);_0!(F);F)

"This is primary data, not necessarily suitable for modified applications"
get_ratea(CaD)= [RateA_Multiplier * CaD["branch"][string(i)]["rate_a"] for i = 1:B]

"""
New wind units are added to bus 8, 65, 116.
Their wind pMax has the ratio 23.90:24.51:72.18 (from the transmission line capacity at that node)
"""
get_windnode() = [8, 65, 116]
get_windPmax() = .54 * [23.9, 24.51, 72.18]

"""
New storage units are added to bus 17, bus 59 (load side), bus 69 and bus 116 (source side)
"""
get_esnode() = [17, 59, 69, 116]
get_esEmax() = [10, 24, 50, 24.] # Assume they (fastest possible) take 2 hours to be fully charged
get_esPmax() = [17., 18., 17., 18.] # based on 1 hour
get_esPmin() = [5., 4., 3., 2.] # based on 1 hour

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
