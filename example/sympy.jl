using ThreeBodyDecay

using SymPy

import SymPy.PyCall
# 
PyCall.pyimport_conda("sympy.physics.wigner", "sympy")
import_from(sympy.physics.wigner)
# wigner_d_small(Sym(1),θ)[1,1] # = cos(θ) d^1_00(θ)

PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
# WignerD

@syms θ
# 
@syms m1 m2 m3 m0
@syms σ1 σ2 σ3
# 
ms = ThreeBodyMasses(m1, m2, m3; m0)
σs = Invariants(ms; σ1, σ2)

# 
cosθ23(σs, ms^2)
cosθ12(σs, ms^2)

# 
@syms z
σ1of2(z, σ2, ms^2)
σ2of1(z, σ1, ms^2)


const J = 3 // 2

function spinparity(p)
    pt = (p[2]..., p[1])
    jpv = str2jp.(pt)
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.SpinTuple,
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.ParityTuple
end

# Ξc(J) -> Ξc(3/2+) [-> Ξc(1/2+) π(0+)] π(0+) 
reaction = "1/2+" => ("1/2+", "0-", "0-")
js, Ps = reaction |> spinparity
tbs = ThreeBodySystem(ms, jv)

Rjp = jp"3/2+"
R(σ) = Sym("R")
# 
dcv = DecayChainsLS(3, R; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
dc = dcv[1, 1]

# 
fullexpression = amplitude(dc[1, 1], σs, (1, 0, 0, 1))

import ThreeBodyDecay: cosθij, cosζ
import ThreeBodyDecay.PartialWaveFunctions: wignerd_doublearg
import ThreeBodyDecay: MandestamTuple, WignerRotation, wr

struct cosHold{T}
    angle::T
end

Sym("θ12") |> typeof

function cosθij(k, σs::MandestamTuple{Sym}, msq)
    i, j, _ = ijk(k)
    θ = Sym(Symbol("θ_", i, j))[1]
    cosHold(θ)
end

label(wr::WignerRotation) = "^$(wr.k)_" *
                            (ispositive(wr) ? "+" : "-") *
                            (iseven(wr) ? "e" : "0")
# 
label(wr(1, 2, 3))
# 
cosζ(wr::WignerRotation{0}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold
cosζ(wr::WignerRotation{2}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold
cosζ(wr::WignerRotation{3}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold

function wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ::cosHold)
    half = 1 / Sym(2)
    WignerD(two_j * half, two_λ1 * half, two_λ2 * half,
        0, cosθ.angle, 0)
end


amplitude(dc, σs, (1, 0, 0, 1); refζs=(1, 2, 3, 1)).doit() |> simplify
# 
sum(itr(js)) do two_λs
    amplitude(dc, σs, two_λs).doit()^2
end |> simplify