using Plots
using ThreeBodyDecay
using LaTeXStrings
pyplot()

#
ms = (Jψ=3.09, p=0.938, Bs=5.366) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.p, m0=ms.Bs;   # masses m1,m2,m3,m0
    two_jps=([1, 1 // 2, 1 // 2, 0] .|> x2,  # twice spin
        ['-', '+', '-', '+'])) # parities
#
# chain 2,3, i.e. (2+3): Λs with the lowest ls, LS
Pc4430_2 = DecayChainLS(2, σ -> BW(σ, 4.330, 0.02); two_s=1 / 2 |> x2, parity='+', tbs=tbs)
Pc4430_3 = DecayChainLS(3, σ -> BW(σ, 4.330, 0.02); two_s=1 / 2 |> x2, parity='-', tbs=tbs)
Pcs = (Pc4430_2, Pc4430_3)
#
# chain-1, i.e. (1+2): Pentaquarks with the lowest ls, LS
ggS = DecayChainLS(1, σ -> 1.0; two_s=1 |> x2, tbs=tbs, parity='-')
ggL = (ggS,)
#
A(σs, two_λs, cs) = sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, (Pcs..., ggL...)))
I(σs, cs) = sum(abs2(A(σs, two_λs, cs)) for two_λs in itr(tbs.two_js))
#
const cs0 = [1.0, -1.0, 20.0 + 0im]
I(σs) = I(σs, cs0) # set the couplings

let
    plot(size=(500, 450),
        xlab=L"m_{p J/\psi}^2\,(\mathrm{GeV})", ylab=L"m_{p J/\psi}^2\,(\mathrm{GeV})")
    #
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ2v = range(tbs.mthsq[2], tbs.sthsq[2], length=100)
    cal = [inphrange((gσ1(σ2, σ3, tbs.ms^2), σ2, σ3), tbs) ? I([gσ1(σ2, σ3, tbs.ms^2), σ2, σ3]) : NaN
           for (σ2, σ3) in Iterators.product(σ2v, σ3v)]
    heatmap!(σ2v, σ3v, cal - cal', colorbar=false, c=:viridis, title="unpolarized intensity")
end

