module WindGen
# TODO add starting time
import Distributions

const m = [ # mean series
    [0.20, 0.22, 0.23, 0.25, 0.28, 0.32, 0.38, 0.45, 0.52, 0.58, 0.62, 0.65,
      0.68, 0.70, 0.73, 0.75, 0.78, 0.80, 0.82, 0.75, 0.60, 0.45, 0.35, 0.28],
    [0.18, 0.20, 0.22, 0.24, 0.29, 0.37, 0.47, 0.58, 0.68, 0.75, 0.80, 0.83,
      0.85, 0.85, 0.83, 0.80, 0.76, 0.70, 0.62, 0.50, 0.42, 0.35, 0.28, 0.22],
    [0.25, 0.27, 0.30, 0.35, 0.45, 0.60, 0.75, 0.82, 0.88, 0.85, 0.78, 0.70,
      0.62, 0.55, 0.50, 0.55, 0.60, 0.68, 0.78, 0.80, 0.70, 0.55, 0.40, 0.30]
]
const Deps = Distributions.MvNormal(
    [1.31 0.6435 0.55
    0.6435 0.79 0.405
    0.55 0.405 1]/200
)
const regression_vec = [0.79, 0.85, 0.87] # inertia coefficients
function m!(W)
    T, N = size(W)
    N === 3 || error("number of wind farms not supported")
    for t = 2:T
        v = rand(Deps)
        for n = 1:N
            k = regression_vec[n]
            v[n] += m[n][t]k + (1-k)W[t-1, n]
        end
        W[t, :] = v
    end
    clamp!(W, 0.0, 1.5)
    W
end

function m24_3()
    W = Matrix{Float64}(undef, 24, 3)
    W[1, :] = [0.187, 0.171, 0.246]
    m!(W)
end

function get_case2383(S)
    Wind = (
        Zone = [1, 2, 3], # Zone Numbers
        N = [ # Node by Zone
            [494, 10, 16, 17, 18],
            [42, 43, 44, 30, 31],
            [105, 64, 1426, 63, 67]
        ],
        PMax = [ # Capacity by Zone
            [2.34, 4.0, 7.2, 10.8, 25.2],
            [4.1, 4.1, 4.1, 8.0, 10.0],
            [4.3, 4.5, 4.95, 6.5, 7.5]
        ],
        S = [ # Sampled trajectories in [0, 1]
            m24_3() for s=1:S
        ] # Wind.S[∀s] has 3 columns. ∀ wind units in Zone 1 uses the 1st column (same reasoning for Zone 2 and 3)
    )
end

end
