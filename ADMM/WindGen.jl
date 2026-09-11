"""
Since I'm doing 2SSP simulation, I need sampler rather than raw data
so I can generate any number of scenes on demand
"""
module WindGen
import Distributions.MvNormal

const (z0, kv, Σ, w1s, w2s, w3s) = ([0.3034779488568444, -0.629419375418374, 0.07677928744799473], [0.9856812720595798, 0.9822456469095393, 0.9663370628399683], [0.027921576704103473 -0.0027188771518694586 0.0010976855707291746; -0.0027188771518694586 0.035576779311823264 -0.0015958186602520468; 0.0010976855707291746 -0.0015958186602520468 0.06444208202855113], (m = -0.8095728876489761, s = 2.9904328434887386), (m = -1.8594780968043005, s = 2.435598925631213), (m = -3.815280149490534, s = 2.6845293048212304));
const ΣD = MvNormal(zeros(3), Σ); # we've got the model param (kv, ΣD)
# Finally, multiply a custom `P_MAX` to convert to real power
function m!(rngW, t, W)
    Taug, N = size(W)
    for t = t+1:Taug
        v = rand(rngW, ΣD)
        @inbounds for n = 1:N
            W[t, n] = kv[n] * W[t-1, n] + v[n]
        end
    end
end
_e(Taug,Ss; NZones=3) = [Matrix{Float64}(undef,Taug,NZones) for _=Ss]
function _S(rngW, S, t; T=16)
    v_of_m = _e(T+1, 1:S)
    Threads.@threads for W = v_of_m
        W[t,:] .= z0 # initialize
        m!(rngW, t, W)
        W[:, 1] .*= w1s.s; W[:, 1] .+= w1s.m
        W[:, 2] .*= w2s.s; W[:, 2] .+= w2s.m
        W[:, 3] .*= w3s.s; W[:, 3] .+= w3s.m
        @. W = logistic(W)
    end
    v_of_m
end
logistic(x) = 1 / (1 + exp(-x))

end
