using Plots
using ThreeBodyDecay
#
# decay Λb ⟶ Jψ p K
ms = (Jψ = 3.09, p=0.938, K = 0.49367, Lb = 5.62) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.K, ms.Lb;   # masses m1,m2,m3,m0
            two_jps=([    1, 1//2,    0,  1//2] .|> x2,  # twice spin
                     [  '-',  '+',  '-',   '+'])) # parities
#
# chain-1, i.e. (2+3): Λs with the lowest ls, LS
Λ1520  = decay_chain(1, (s,σ)->BW(σ, 1.5195, 0.017); two_s = 3/2|>x2, tbs=tbs)
Λ1690  = decay_chain(1, (s,σ)->BW(σ, 1.685,  0.050); two_s = 1/2|>x2, tbs=tbs)
Λ1810  = decay_chain(1, (s,σ)->BW(σ, 1.80,   0.090); two_s = 5/2|>x2, tbs=tbs)
Λs = (Λ1520,Λ1690,Λ1810)
#
# chain-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
Pc4312 = decay_chain(3, (s,σ)->BW(σ, 4.312, 0.015); two_s = 1/2|>x2, tbs=tbs)
Pc4440 = decay_chain(3, (s,σ)->BW(σ, 4.440, 0.010); two_s = 1/2|>x2, tbs=tbs)
Pc4457 = decay_chain(3, (s,σ)->BW(σ, 4.457, 0.020); two_s = 3/2|>x2, tbs=tbs)
Pcs = (Pc4312,Pc4440,Pc4457)
#
A(σs,two_λs,cs) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(cs, (Λs...,Pcs...)))
I(σs,cs) = sum(abs2(A(σs,two_λs,cs)) for two_λs in itr(tbs.two_js))
#
I(σs) = I(σs,[1, 1.1, 3.9im, 2.2, 2.1im, -0.3im]) # set the couplings
#
dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
@show I(dpp.σs)

function plot_something(what; method=heatmap)
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ1v = range(tbs.mthsq[1], tbs.sthsq[1], length=100)
    cal = [Kibble31(σ3,σ1,tbs.msq) < 0.0 ? what([σ1,gσ2(σ3,σ1,tbs.msq),σ3]) : NaN
        for (σ3,σ1) in Iterators.product(σ3v,σ1v)]
    method(σ1v, σ3v, cal, c=:lime_grad, colorbar=false)
end

@time plot_something(I)

using Plots
using LaTeXStrings
pyplot()

let
    plot(size=(1000,350), layout=grid(1,2))
    plot!(sp=1, border31(tbs),
        xlab=L"\sigma_1\,(\mathrm{GeV})", ylab=L"\sigma_3\,(\mathrm{GeV})",
        lab="", l=(2,:black))
    #
    plot!(sp=2, border12(tbs),
        xlab=L"\sigma_2\,(\mathrm{GeV})", ylab=L"\sigma_1\,(\mathrm{GeV})",
        lab="", l=(2,:black))
end
savefig(joinpath("example","plot","border31_12.pdf"))
savefig(joinpath("example","plot","border31_12.png"))

let
    # generate data
    plot(size=(1000,350), layout=grid(1,2),
        xlab=L"\sigma_1\,(\mathrm{GeV})", ylab=L"\sigma_3\,(\mathrm{GeV})")
    σ1v,σ2v,σ3v = flatDalitzPlotSample(tbs; Nev = 10000)
    scatter!(sp=1, σ1v,σ3v, lab="", markerstrokewidth=0.5)
    # weight with amplitude
    weights = [I(σs) for σs in zip(σ1v,σ2v,σ3v)]
    # weighted histogram
    histogram2d!(sp=2, σ1v,σ3v, weights=weights, bins=60,
        c=:viridis)
end
savefig(joinpath("example","plot","dalitz31.pdf"))
savefig(joinpath("example","plot","dalitz31.png"))
