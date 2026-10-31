"""
random wind generation
"""
module WindGen
import Random, Distributions.MvNormal
const T = 12 # num of look-ahead periods
const #= fitted AR model =# (kv, Σ) = ([0.9793411235799021, 0.9820667500156772, 0.9774222552049533], [0.004090916685855472 0.0007655455638262616 0.00017323193031670899; 0.0007655455638262616 0.003455090260169231 0.00012804507593629783; 0.00017323193031670899 0.00012804507593629783 0.004461276600803643])
const #= z_score =# z0_default = [-1.0204188052610943, -0.7449981411563563, -0.32278370597619005]
const Wμσvec = [(m = -0.1940367284110205, s = 2.459001343589591), (m = -1.1467751498561531, s = 3.031566738500728), (m = -2.430637986045454, s = 2.3876609351738525)]
const zmat = zeros(1+T, length(kv)) # Notice that Julia is 1-based index
const ΣD = MvNormal(Σ)

logistic(x) = 1 / (1 + exp(-x))

function step!(zmat, t, rng)
    i = t+1 # julia index
    v = rand(rng, ΣD)
    v .+= (kv .* zmat[i-1, :])
    zmat[i, :] .= v
end
reveal!(zmat, t, rng) = if t == 0
    zmat[1, :] .= z0_default
else
    step!(zmat, t, rng)
end
revealall!(zmat, t, rng) = for t = t:size(zmat, 1)-1
    step!(zmat, t, rng)
end


function zmat2Wmat(zmat)
    Wmat = similar(zmat)
    for (z,nt) = enumerate(Wμσvec), i = 1:size(zmat, 1)
        Wmat[i, z] = logistic(zmat[i, z] * nt.s + nt.m)
    end
    Wmat
end
fill_Wscnvec_once!(Wscnvec, zmat, t, S) = for s=1:S # this is an initial Wscnvec
    revealall!(zmat, t, Random.Xoshiro(767234t+s))
    Wscnvec[s] = zmat2Wmat(zmat)
end

end


