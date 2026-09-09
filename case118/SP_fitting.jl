# This file is to demonstrate how to simulate wind power series from wind farms in multi-locations
# For theories, refer to wind_farms_model.png

using GLMakie
# Makie.update_theme!(
    # fonts = (
        # regular = "Microsoft YaHei",
        # bold = "Microsoft YaHei Bold"
    # )
# )

import Statistics, LinearAlgebra, Distributions
lin_coeff(v) = sum(v[i]v[i-1] for i=2:length(v)) / sum(v[i-1]^2 for i=2:length(v))
logit(x) = log(x / (1 - x));
logistic(x) = 1 / (1 + exp(-x));
# load the source time series data, e.g. from Inner-Mongolia wind farms (at 3 different locations)
w1 = [25.18, 25.6, 25.72, 26.94, 27.67, 26.98, 20.99, 15.7, 12.19, 10.35, 7.02, 6.29, 4.85, 4.7, 5.63, 6.7, 7.17, 7.74, 6.84, 7.85, 7.03, 6.64, 6.81, 7.71, 9.15, 11.08, 9.18, 8.89, 9.9, 7.74, 5.22, 3.17, 3.29, 6.32, 7.91, 8.32, 9.14, 7.35, 6.03, 3.3, 2.45, 2.63, 2.06, 2.47, 4.57, 4.83, 4.61, 5.75, 5.35, 4.78, 5.25, 4.01, 4.38, 5.44, 4.29, 4.2, 4.29, 5.18, 5.05, 4.58, 2.67, 2.76, 3.85, 3.62, 2.69, 2.72, 1.11, 0.0, 0.0, 0.0, 0.52, 0.62, 0.44, 0.23, 0.17, 0.29, 0.22, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 1.05, 1.13, 2.03, 8.35];
w2 = [2.331, 2.405, 1.592, 1.909, 2.391, 2.856, 3.47, 4.142, 5.997, 4.733, 2.896, 2.599, 2.831, 3.445, 4.76, 5.669, 5.374, 6.053, 5.358, 4.416, 4.771, 5.251, 5.253, 7.339, 9.25, 11.121, 13.807, 14.083, 12.612, 11.719, 11.505, 15.048, 15.475, 15.658, 14.403, 13.28, 11.692, 10.525, 8.83, 8.395, 8.102, 7.687, 7.301, 7.057, 6.615, 8.591, 9.864, 8.417, 7.153, 6.919, 6.638, 5.202, 4.947, 5.191, 4.626, 3.916, 4.088, 4.025, 5.454, 6.682, 8.37, 9.797, 10.576, 7.634, 5.584, 5.432, 5.622, 5.986, 6.14, 5.637, 6.126, 5.073, 3.74, 3.952, 3.757, 2.389, 1.634, 0.871, 0.123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
w3 = [0.654, 0.717, 0.694, 0.927, 0.533, 0.493, 0.513, 0.129, 0.002, 0.145, 0.042, 0.0, 0.04, 0.589, 1.043, 1.282, 1.748, 1.666, 2.279, 3.746, 7.701, 7.475, 6.457, 7.332, 7.955, 8.435, 8.107, 8.48, 7.636, 7.174, 6.475, 4.927, 3.983, 3.88, 3.481, 2.88, 2.492, 1.862, 0.958, 0.886, 0.663, 0.404, 0.373, 0.404, 1.949, 3.139, 3.032, 3.411, 2.822, 3.646, 3.014, 2.494, 1.951, 2.103, 1.904, 1.279, 0.743, 0.263, 0.0, 0.0, 0.054, 0.377, 0.125, 0.002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.013, 1.583, 3.88, 5.544, 7.013, 7.788, 8.426, 10.773, 11.949, 13.418, 10.083, 7.977, 8.227, 6.662, 5.063, 5.316, 5.562, 5.702, 4.727, 4.983, 7.301, 16.446, 16.738, 14.175, 10.022, 8.821];
for k=(w1,w2,w3) k ./= maximum(k) end # result in [0,1] values
for k=(w1,w2,w3), (i,e)=enumerate(k) k[i]=clamp(e, 0.001, 0.999) end # result in (0,1) values
for k=(w1,w2,w3) @. k = logit(k) end # result in values in the \mathbb{R} axis
w1s = (m = Statistics.mean(w1), s = Statistics.std(w1; corrected=false));
w2s = (m = Statistics.mean(w2), s = Statistics.std(w2; corrected=false));
w3s = (m = Statistics.mean(w3), s = Statistics.std(w3; corrected=false));
z1 = (w1 .- w1s.m) ./ w1s.s;
z2 = (w2 .- w2s.m) ./ w2s.s;
z3 = (w3 .- w3s.m) ./ w3s.s; # Here we have derived z-score
z0 = [z1[1], z2[1], z3[1]]; # initial-time vector for the 3 wind farms
k1 = lin_coeff(z1);
k2 = lin_coeff(z2);
k3 = lin_coeff(z3); # slopes
kv = [k1, k2, k3]; # collect the slopes into a vector
Σ, mat = zeros(3,3), zeros(3,3)
for i = 2:length(z1)
    v = [z1[i] - k1 * z1[i-1], z2[i] - k2 * z2[i-1], z3[i] - k3 * z3[i-1]]
    LinearAlgebra.mul!(mat, v, v')
    Σ .+= mat
end
Σ ./= length(2:length(z1));
ΣD = Distributions.MvNormal(zeros(3), Σ); # we've got the model param (kv, ΣD)
# Here simulation begins
SimL = 12
zrnd = zeros(SimL,3);
zrnd[1, :] .= z0 # initial-time
for i = 2:SimL
    tmpv = rand(ΣD)
    tmpv .+= (kv .* zrnd[i-1, :])
    zrnd[i, :] .= tmpv
end
zrnd[:, 1] .*= w1s.s;
zrnd[:, 1] .+= w1s.m;
zrnd[:, 2] .*= w2s.s;
zrnd[:, 2] .+= w2s.m;
zrnd[:, 3] .*= w3s.s;
zrnd[:, 3] .+= w3s.m;
@. zrnd = logistic(zrnd);
# Finally, we still need 1 more step: `*P_MAX` to convert to real power

f = Figure();
ax1 = Axis(f[1, 1])
xtks = 1:SimL
lines!(ax1, xtks, zrnd[:,1], color = :red)
lines!(ax1, xtks, zrnd[:,2], color = :blue)
lines!(ax1, xtks, zrnd[:,3], color = :green)
f
