"""
Since I'm doing 2SSP simulation, I need sampler rather than raw data so I can generate any number of scenes on demand
"""
module WindGen
import Distributions.MvNormal

# mean series, starting from 16:00, each row is a wind farm, 3*(1+16); where (15min)*16, The 1st column is already-revealed
# if you want to do 5min precision, do interpolation
const mu_06 = Matrix(transpose([
    0.642  0.668  0.704  0.731  0.747  0.739  0.712  0.681  0.659  0.645  0.632  0.616  0.592  0.561  0.538  0.519  0.504
    0.703  0.685  0.662  0.648  0.655  0.679  0.717  0.751  0.773  0.781  0.769  0.744  0.718  0.691  0.663  0.639  0.621
    0.781  0.803  0.827  0.846  0.859  0.872  0.885  0.893  0.887  0.872  0.848  0.821  0.794  0.768  0.743  0.721  0.703
]));
const Deps = MvNormal(
    [1.31 0.6435 0.55
    0.6435 0.79 0.405
    0.55 0.405 1]/200/8 # The last coeff is for smaller variance for intra-day
);
const regression_vec = 0.12 .+ [0.79, 0.85, 0.87] # inertia coefficients

function m!(WmuMat, t, W)
    Taug, N = size(W)
    N === 3 || error("number of wind farms not supported")
    @inbounds for t = t+1:Taug
        v = rand(Deps)
        for n = 1:N
            k = regression_vec[n]
            v[n] += WmuMat[t, n]k + (1-k)W[t-1, n]
        end
        W[t, :] = v
    end
    clamp!(W, 0.0, 1.5)
end

_e(Taug,Ss; NZones=3) = [Matrix{Float64}(undef,Taug,NZones) for _=Ss]
function _S(S, t; T=16, WmuMat=mu_06)
    v_of_m = _e(T+1, 1:S)
    Threads.@threads for W = v_of_m
        W[t,:] .= WmuMat[t,:] # initialize
        m!(WmuMat, t, W)
    end
    v_of_m
end

end
