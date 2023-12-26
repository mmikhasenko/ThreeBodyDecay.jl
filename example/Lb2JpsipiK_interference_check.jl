using Plots
using Plots.PlotMeasures: mm
using LaTeXStrings
#
using ThreeBodyDecay
#
# decay Λb ⟶ Jψ p K
const mJψ = 3.09
const mp = 0.938
const mK = 0.49367
const mLb = 5.62

# create two-body system
tbs = ThreeBodySystem(mJψ, mp, mK; m0=mLb, # masses m1,m2,m3,m0
    two_js=ThreeBodySpins(2, 1, 0; two_h0=1))
#
Ps = ['-', '+', '-', '+']
# chain-1, i.e. (2+3): Λs with the lowest ls, LS
Λ1520 = DecayChainLS(1, σ -> BW(σ, 1.5195, 0.017); two_s=3 / 2 |> x2, tbs, Ps)
Λ1690 = DecayChainLS(1, σ -> BW(σ, 1.685, 0.050); two_s=1 / 2 |> x2, tbs, Ps)
Λ1810 = DecayChainLS(1, σ -> BW(σ, 1.80, 0.090); two_s=5 / 2 |> x2, tbs, Ps)
Λs = (Λ1520, Λ1690, Λ1810)
#
# chain-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
Pc4312 = DecayChainLS(3, σ -> BW(σ, 4.312, 0.015); two_s=1 / 2 |> x2, tbs, Ps)
Pc4440 = DecayChainLS(3, σ -> BW(σ, 4.440, 0.010); two_s=1 / 2 |> x2, tbs, Ps)
Pc4457 = DecayChainLS(3, σ -> BW(σ, 4.457, 0.020); two_s=3 / 2 |> x2, tbs, Ps)
Pcs = (Pc4312, Pc4440, Pc4457)
#
A(σs, two_λs, cs) = sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, (Λs..., Pcs...)))
I(σs, cs) = sum(abs2(A(σs, two_λs, cs)) for two_λs in itr(tbs.two_js))
#
Ini(σs, cs) = sum(abs2(sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, Λs))) for two_λs in itr(tbs.two_js)) +
              sum(abs2(sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, Pcs))) for two_λs in itr(tbs.two_js))
#
coupligns = [1, 1.1, 3.9im, 2.2, 2.1im, -0.3im]
I(σs) = I(σs, coupligns) # set the couplings
Ini(σs) = Ini(σs, coupligns) # set the couplings
#
dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
@show I(dpp.σs)
@show Ini(dpp.σs)

# generate data
data = let
    #
    σsv = flatDalitzPlotSample(tbs.ms; Nev=10000)
    cosθ12v = cosθ12.(σsv, Ref(tbs.ms^2))
    #
    (σsv, cosθ12v,
        weights=I.(σsv), weights_ni=Ini.(σsv))
end

;
let bins = 60, seriestype = :stephist
    plot(size=(900, 350), layout=grid(1, 2),
        xlab=L"\cos\,\theta_{12}", bottom_margin=5mm)
    plot!(sp=1, data.cosθ12v; lab="ph.ps.", seriestype, bins)
    plot!(sp=2, data.cosθ12v; weights=data.weights, lab="|L+P|^2", seriestype, bins)
    plot!(sp=2, data.cosθ12v; weights=data.weights_ni, lab="|L|^2+|P|^2", seriestype, bins)
end