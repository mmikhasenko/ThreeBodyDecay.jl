using ThreeBodyDecay
using LaTeXStrings
using Plots

theme(:wong, frame=:box, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right)

const mD = 1.86483
const mπ0 = 0.13498
const mDstar = 2.00685 # 2 GeV
const ΓDstar = 30e-6 # 30 keV
ms = ThreeBodyMasses(mD, mπ0, mD; m0=mD + mDstar + 0.01)

# amplitude

BWD(σ) = BW(σ, mDstar, 100 * ΓDstar)
A(σs) = BWD(σs.σ1) + BWD(σs.σ3)
A(σs, ν) = BWD(σs.σ1) * wignerd(1, ν, 0, cosθ23(σs, ms^2)) -
           BWD(σs.σ3) * sum(wignerd(1, ν, λ, cosζ31_for0(σs, ms^2)) *
                            wignerd(1, λ, 0, cosθ12(σs, ms^2)) for λ = -1:1)

# intensity
I(σs) = sum(abs2(A(σs, ν)) for ν = -1:1)

let
    plot(ms, I; iσx=1, iσy=3)
    plot!(xlab=L"m(D)\,\,[\mathrm{GeV}]", ylab=L"m(Dπ)\,\,[\mathrm{GeV}]")
end