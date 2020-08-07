using Plots
using ThreeBodyDecay
using LaTeXStrings
pyplot()

#
ms = (Jψ = 3.09, p=0.938, Bs = 5.366) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.p, ms.Jψ, ms.p,  ms.Bs;   # masses m1,m2,m3,m0
            two_jps=([  1//2,   1, 1//2,    0] .|> x2,  # twice spin
                     [   '+', '-',  '-',  '+'])) # parities
#
# chain 2,3, i.e. (2+3): Λs with the lowest ls, LS
Pc4430_1  = decay_chain(1, (s,σ)->BW(σ, 4.330, 0.02); two_s = 1/2|>x2, parity='+', tbs=tbs)
Pc4430_3  = decay_chain(3, (s,σ)->BW(σ, 4.330, 0.02); two_s = 1/2|>x2, parity='-', tbs=tbs)
Pcs = (Pc4430_1, Pc4430_3)
#
# chain-1, i.e. (1+2): Pentaquarks with the lowest ls, LS
ggS = decay_chain(2, (s,σ)->1.0; two_s = 1|>x2, tbs=tbs, parity='-')
ggL = (ggS,)
#
A(σs,two_λs,cs) = sum(c*amplitude(σs,two_λs,dc) for (c, dc) in zip(cs, (Pcs..., ggL...)))
I(σs,cs) = sum(abs2(A(σs,two_λs,cs)) for two_λs in itr(tbs.two_js))
#
const cs0 = [1.0, -1.0, 20.0+0im]
I(σs) = I(σs,cs0) # set the couplings
#
# dpp = randomPoint(tbs) # just a random point of the Dalitz Plot
# @show I(dpp.σs)
#
let
    plot(size=(500,450),
        xlab=L"m_{p J/\psi}^2\,(\mathrm{GeV})", ylab=L"m_{p J/\psi}^2\,(\mathrm{GeV})")
    #
    σ3v = range(tbs.mthsq[3], tbs.sthsq[3], length=100)
    σ1v = range(tbs.mthsq[1], tbs.sthsq[1], length=100)
    cal = [Kibble31(σ3,σ1,tbs.msq) < 0.0 ? I([σ1,gσ2(σ3,σ1,tbs.msq),σ3]) : NaN
        for (σ3, σ1) in Iterators.product(σ3v,σ1v)]
    heatmap!(σ1v, σ3v, cal - cal', colorbar=false, c=:viridis, title="unpolarized intensity")
end
