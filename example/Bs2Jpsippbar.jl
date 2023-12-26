using Plots
using ThreeBodyDecay
using LaTeXStrings

#
ms = (Jψ=3.09, p=0.938, Bs=5.366) # masses of the particles
# create two-body system
tbs = ThreeBodySystem(ms.Jψ, ms.p, ms.p; m0=ms.Bs,   # masses m1,m2,m3,m0
    two_js=ThreeBodySpins(2, 1, 1; two_h0=0)  # twice spin
)

Ps = ThreeBodyParities('+', '-', '+'; P0='+') # parities of the particles
#
# chain 2,3, i.e. (2+3): Λs with the lowest ls, LS
Pc4430_2 = DecayChainLS(2, σ -> BW(σ, 4.330, 0.02); two_s=1 / 2 |> x2, parity='+', tbs=tbs, Ps)
Pc4430_3 = DecayChainLS(3, σ -> BW(σ, 4.330, 0.02); two_s=1 / 2 |> x2, parity='-', tbs=tbs, Ps)
Pcs = (Pc4430_2, Pc4430_3)
#
# chain-1, i.e. (1+2): Pentaquarks with the lowest ls, LS
ggS = DecayChainLS(1, σ -> 1.0; two_s=1 |> x2, tbs=tbs, parity='-', Ps)
ggL = (ggS,)
#
A(σs, two_λs, cs) = sum(c * amplitude(dc, σs, two_λs) for (c, dc) in zip(cs, (Pcs..., ggL...)))
I(σs, cs) = sum(abs2(A(σs, two_λs, cs)) for two_λs in itr(tbs.two_js))
#
const cs0 = [1.0, -1.0, 20.0 + 0im]
I(σs) = I(σs, cs0) # set the couplings

let
    plot(size=(500, 450))
    plot!(tbs.ms, I; iσx=2, iσy=3)
    plot!(border23(tbs.ms), lc=2, lw=5, lab="")
    plot!(
        xlab=L"m_{p J/\psi}^2\,(\mathrm{GeV})",
        ylab=L"m_{p J/\psi}^2\,(\mathrm{GeV})")
end

