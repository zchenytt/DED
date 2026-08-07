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

"energy storage"
const EEini = [5.9, 14.16, 29.5, 9.84]

const B, N, L = 186, 118, 99-length(EVnode)
const Load_Multiplier = 1.0 # to scale up the load to match the new wind capacity
const RateA_Multiplier = 1.0
const LookAheadHours = 4 # We always look ahead 4 hours; Then we have 2 possibilities: (5min)48 _or_ (15min)16
const MyRateA = [0.03, 0.49, 1.29, 0.8, 1.13, 0.73, 1.32, 4.29, 1.32, 0.94, 1.09, 0.84, 0.23, 0.05, 0.56, 0.74, 0.6, 0.43, 0.52, 0.5, 1.1, 0.37, 1.01, 0.45, 0.3, 0.59, 0.14, 0.13, 0.19, 0.6, 1.04, 1.4, 0.81, 0.13, 0.13, 2.4, 4.42, 0.78, 0.78, 0.35, 0.43, 0.12, 0.16, 0.48, 0.35, 0.09, 0.31, 0.31, 0.26, 0.84, 3.31, 1.3, 1.16, 2.07, 1.04, 0.92, 0.68, 0.59, 0.26, 0.4, 0.12, 0.21, 0.31, 0.1, 0.28, 0.77, 0.77, 0.28, 0.24, 0.72, 0.89, 0.4, 0.24, 0.04, 0.66, 0.64, 0.08, 0.64, 0.26, 0.49, 0.58, 0.31, 0.38, 0.44, 0.36, 0.38, 0.45, 1.15, 1.25, 1.53, 0.41, 0.19, 4.35, 4.35, 2.5, 1.91, 6.73, 2.19, 2.19, 0.55, 0.45, 4.7, 0.68, 4.02, 0.84, 0.77, 5.62, 1.28, 0.22, 0.11, 0.27, 0.16, 0.06, 0.46, 0.38, 1.65, 0.2, 0.06, 2.9, 0.45, 0.8, 0.15, 0.26, 0.13, 0.33, 8.51, 8.51, 2.35, 1.67, 0.62, 0.87, 0.52, 0.19, 0.08, 0.62, 0.44, 0.18, 0.97, 1.86, 0.93, 0.67, 0.21, 1.0, 0.75, 0.79, 0.84, 0.85, 1.29, 0.29, 1.2, 1.37, 1.47, 1.37, 0.46, 0.45, 1.23, 1.24, 1.18, 1.08, 0.58, 0.38, 0.41, 1.33, 0.56, 0.27, 0.35, 0.56, 0.47, 0.18, 0.28, 0.31, 0.23, 0.29, 0.6, 0.22, 0.63, 0.63, 0.59, 0.54, 0.22, 0.16, 0.14, 17.95, 0.19, 0.97, 0.67]

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
get_windPmax() = .457 * [23.9, 24.51, 72.18]

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
function get_myGPmax(CaD, Ggvec)
    GPmax = get_GPmax(CaD, Ggvec)
    for (i, P) in enumerate(GPmax)
        P > 1 && (GPmax[i] /= 3)
    end
    GPmax
end
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
