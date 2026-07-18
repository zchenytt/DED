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
const Load_Multiplier = 1.94 # to scale up the load to match the new wind capacity
const RateA_Multiplier = 1.0

# Nominal load 5min * (1+36)
const l00 = ( # 00:00 - 03:00 load types
    [0.35, 0.34, 0.33, 0.32, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.20, 0.20, 0.19, 0.19, 0.18, 0.18, 0.17, 0.17, 0.17, 0.16, 0.16, 0.16, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.16, 0.16, 0.16, 0.16],
    [0.13, 0.13, 0.13, 0.12, 0.11, 0.11, 0.11, 0.12, 0.13, 0.13, 0.13, 0.12, 0.11, 0.11, 0.11, 0.12, 0.13, 0.13, 0.13, 0.12, 0.11, 0.11, 0.11, 0.12, 0.13, 0.13, 0.13, 0.12, 0.11, 0.11, 0.11, 0.12, 0.13, 0.13, 0.13, 0.12, 0.11],
    [0.85, 0.85, 0.84, 0.84, 0.85, 0.85, 0.86, 0.86, 0.85, 0.85, 0.84, 0.84, 0.84, 0.85, 0.85, 0.85, 0.84, 0.84, 0.85, 0.85, 0.86, 0.86, 0.85, 0.85, 0.84, 0.84, 0.84, 0.85, 0.85, 0.85, 0.84, 0.84, 0.85, 0.85, 0.86, 0.86, 0.85]
);
const l06 = ( # 06:00 - 09:00 load types
    [0.26, 0.26, 0.26, 0.28, 0.31, 0.34, 0.38, 0.41, 0.44, 0.45, 0.46, 0.48, 0.51, 0.55, 0.58, 0.61, 0.63, 0.64, 0.65, 0.66, 0.68, 0.70, 0.73, 0.76, 0.78, 0.79, 0.79, 0.79, 0.80, 0.81, 0.83, 0.84, 0.84, 0.83, 0.81, 0.79, 0.76],
    [0.15, 0.16, 0.18, 0.20, 0.23, 0.27, 0.32, 0.38, 0.44, 0.50, 0.55, 0.60, 0.65, 0.70, 0.74, 0.78, 0.81, 0.84, 0.86, 0.88, 0.89, 0.90, 0.91, 0.92, 0.92, 0.93, 0.93, 0.94, 0.94, 0.94, 0.95, 0.95, 0.95, 0.95, 0.96, 0.96, 0.96],
    [0.85, 0.85, 0.84, 0.83, 0.82, 0.80, 0.78, 0.77, 0.76, 0.77, 0.79, 0.82, 0.84, 0.86, 0.88, 0.89, 0.90, 0.91, 0.91, 0.92, 0.92, 0.91, 0.91, 0.92, 0.92, 0.92, 0.93, 0.93, 0.93, 0.93, 0.94, 0.94, 0.94, 0.94, 0.94, 0.95, 0.95]
);
const l15 = ( # 15:00 - 18:00 load types
    [0.55, 0.56, 0.57, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.73, 0.75, 0.78, 0.80, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.99, 1.00, 1.00, 0.99, 0.99, 0.98, 0.98, 0.97, 0.97, 0.96, 0.96],
    [0.95, 0.95, 0.96, 0.96, 0.95, 0.95, 0.94, 0.94, 0.94, 0.93, 0.93, 0.92, 0.92, 0.91, 0.90, 0.89, 0.88, 0.86, 0.84, 0.82, 0.80, 0.77, 0.75, 0.72, 0.69, 0.67, 0.65, 0.63, 0.61, 0.59, 0.57, 0.55, 0.54, 0.53, 0.52, 0.51, 0.50],
    [0.94, 0.94, 0.95, 0.95, 0.94, 0.94, 0.95, 0.95, 0.94, 0.94, 0.95, 0.95, 0.94, 0.94, 0.93, 0.93, 0.92, 0.91, 0.90, 0.91, 0.92, 0.93, 0.94, 0.94, 0.95, 0.95, 0.94, 0.94, 0.94, 0.95, 0.95, 0.94, 0.94, 0.95, 0.95, 0.94, 0.94]
);
get_load_type() = rand(eachindex(l00), L)
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
get_esEmax() = [34., 35.] # a unit length of period default to 15min
get_esPmax() =  [8.,  9.]

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
