using Plots
using LaTeXStrings
#
using ThreeBodyDecay
#
# decay Λb ⟶ Jψ p K
ms = (Jψ=3.09, p=0.938, K=0.49367, Lb=5.62) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.K; m0=ms.Lb,   # masses m1,m2,m3,m0
    two_jps=([1, 1 // 2, 0, 1 // 2] .|> x2,  # twice spin
        ['-', '+', '-', '+'])) # parities
#
# chain-1, i.e. (2+3): Λs with the lowest ls, LS
Λ1520 = DecayChainLS(1, σ -> BW(σ, 1.5195, 0.017); two_s=3 / 2 |> x2, tbs=tbs)
Λ1690 = DecayChainLS(1, σ -> BW(σ, 1.685, 0.050); two_s=1 / 2 |> x2, tbs=tbs)
Λ1810 = DecayChainLS(1, σ -> BW(σ, 1.80, 0.090); two_s=5 / 2 |> x2, tbs=tbs)
Λs = (Λ1520, Λ1690, Λ1810)
#
# chain-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
Pc4312 = DecayChainLS(3, σ -> BW(σ, 4.312, 0.015); two_s=1 / 2 |> x2, tbs=tbs)
Pc4440 = DecayChainLS(3, σ -> BW(σ, 4.440, 0.010); two_s=1 / 2 |> x2, tbs=tbs)
Pc4457 = DecayChainLS(3, σ -> BW(σ, 4.457, 0.020); two_s=3 / 2 |> x2, tbs=tbs)
Pcs = (Pc4312, Pc4440, Pc4457)
#
A(σs, two_λs, cs) = sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, (Λs..., Pcs...)))
I(σs, cs) = sum(abs2(A(σs, two_λs, cs)) for two_λs in itr(tbs.two_js))
#
Ini(σs, cs) = sum(abs2(sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, Λs))) for two_λs in itr(tbs.two_js)) +
              sum(abs2(sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, Pcs))) for two_λs in itr(tbs.two_js))
#
I(σs) = I(σs, [1, 1.1, 3.9im, 2.2, 2.1im, -0.3im]) # set the couplings
Ini(σs) = Ini(σs, [1, 1.1, 3.9im, 2.2, 2.1im, -0.3im]) # set the couplings
#
dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
@show I(dpp.σs)
@show Ini(dpp.σs)

let
    # generate data
    plot(size=(900, 350), layout=grid(1, 2),
        xlab=L"\cos\,\theta_{12}")
    #
    σ1v, σ2v, σ3v = flatDalitzPlotSample(tbs; Nev=10000)
    cosθ12v = [cosθ12(σs, tbs.msq) for σs in zip(σ1v, σ2v, σ3v)]
    #
    weights = [I(σs) for σs in zip(σ1v, σ2v, σ3v)]
    weights_ni = [Ini(σs) for σs in zip(σ1v, σ2v, σ3v)]
    #
    plot!(sp=1, cosθ12v, bins=60, seriestype=:stephist, lab="ph.ps.")
    plot!(sp=2, cosθ12v, weights=weights, bins=60, seriestype=:stephist, lab="|L+P|^2")
    plot!(sp=2, cosθ12v, weights=weights_ni, bins=60, seriestype=:stephist, lab="|L|^2+|P|^2")
end
