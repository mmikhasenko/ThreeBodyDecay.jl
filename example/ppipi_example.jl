using Plots
using ThreeBodyDecay
#
# decay Λb ⟶ Jψ p K
ms = (π = 0.14, p=0.938, W = 2.4) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.p, ms.π, ms.π, ms.W;   # masses m1,m2,m3,m0
            two_jps=([ 1//2,    0,    0,  3//2] .|> x2,  # twice spin
                     [  '+',  '-',  '-',   '+'])) # parities
#
# chain-1, i.e. (2+3): Λs with the lowest ls, LS
Δ1232  = [decay_chain(ch, (s,σ)->BW(σ, 1.232, 0.13); two_s = 3/2|>x2, tbs=tbs) for ch=[2,3]]
Δs = (Δ1232...,)
#
# chain-3, i.e. (1+2): Pentaquarks with the lowest ls, LS
ρ770 = decay_chain(1, (s,σ)->BW(σ, 0.775, 0.15); two_s = 1|>x2, tbs=tbs, parity='-')
f2 = decay_chain(1, (s,σ)->BW(σ, 1.26, 0.186); two_s = 2|>x2, tbs=tbs)
ππL = (ρ770,f2)
#
A(σs,two_λs,cs) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(cs, (Δs..., ππL...)))
I(σs,cs) = sum(abs2(A(σs,two_λs,cs)) for two_λs in itr(tbs.two_js))
#
I(σs) = I(σs,[1, 1, 0.9, 1.1im]) # set the couplings
#
dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
@show I(dpp.σs)

using LaTeXStrings
pyplot()

let
    plot(size=(500,450),
        xlab=L"m_{p\pi^+}^2\,(\mathrm{GeV})", ylab=L"m_{p\pi^-}^2\,(\mathrm{GeV})")
    #
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ2v = range(tbs.mthsq[2], tbs.sthsq[2], length=100)
    cal = [Kibble23(σ2,σ3,tbs.msq) < 0.0 ? I([gσ1(σ2,σ3,tbs.msq),σ2,σ3]) : NaN
        for (σ2,σ3) in Iterators.product(σ2v,σ3v)]
    heatmap!(σ2v, σ3v, cal, colorbar=false, c=:viridis, title="unpolarized intensity")
end

let
    # generate data
    plot(size=(500,450),
        xlab=L"m_{p\pi^+}^2\,(\mathrm{GeV})", ylab=L"m_{p\pi^-}^2\,(\mathrm{GeV})")
    σ1v,σ2v,σ3v = flatDalitzPlotSample(tbs; Nev = 10000)
    # weight with amplitude
    weights = [I(σs) for σs in zip(σ1v,σ2v,σ3v)]
    # weighted histogram
    histogram2d!(σ3v, σ2v, weights=weights, bins=60,
        c=:viridis, title="weighted MC")
end
