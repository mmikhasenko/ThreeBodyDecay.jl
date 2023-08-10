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


collect("1/2+" => ("1/2+", "0-", "0-"))
# Ξc(J) -> Ξc(3/2+) [-> Ξc(1/2+) π(0+)] π(0+) 
reaction = "1/2+" => ("1/2+", "0-", "0-")
js, Ps = reaction |> spinparity
tbs = ThreeBodySystem(ms, js)

Rjp = jp"3/2+"
R(σ) = Sym("R")
# 
dc = DecayChainLS(3, R; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)

# 
fullexpression = amplitude(dc, σs, (1, 0, 0, 1))

import ThreeBodyDecay: cosθij, cosζ
import ThreeBodyDecay.PartialWaveFunctions: wignerd_doublearg
import ThreeBodyDecay: MassTuple, MandestamTuple, WignerRotation, wr

struct cosHold{T,X}
    angle::T
    expession::X
end

struct StickySymTuple{X,N}
    data::NamedTuple{X,NTuple{N,Sym}}
end
Base.getproperty(mnt::StickySymTuple, name::Symbol) = getproperty(getfield(mnt, :data), name)
Base.getindex(mnt::StickySymTuple, i::Int) = getfield(mnt, :data)[i]

sσs = StickySymTuple(σs)

function cosθij(k, σs::StickySymTuple{(:σ1, :σ2, :σ3),3}, msq)
    i, j, _ = ijk(k)
    θ = Sym(Symbol("θ_", i, j))[1]
    cosHold(θ, cosθij(k, getfield(σs, :data), msq))
end

cosθ12(σs, ms^2)
cosθ12(σs |> StickySymTuple, ms^2)

label(wr::WignerRotation) = "^$(wr.k)_" *
                            (ispositive(wr) ? "+" : "-") *
                            (iseven(wr) ? "e" : "0")
# 
for N in (0, 2, 3)
    eval(:(
        function cosζ(wr::WignerRotation{$(N)},
            σs::StickySymTuple{(:σ1, :σ2, :σ3),3}, msq)
            ζ = Sym("ζ" * label(wr))
            return cosHold(ζ, cosζ(wr, getfield(σs, :data), msq))
        end
    ))
end

cosζ(wr(1, 2, 1), σs, ms^2)
cosζ(wr(1, 2, 1), σs |> StickySymTuple, ms^2)

function wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ::cosHold)
    half = 1 / Sym(2)
    WignerD(two_j * half, two_λ1 * half, two_λ2 * half,
        0, cosθ.angle, 0)
end

amplitude(dc, σs |> StickySymTuple, (1, 0, 0, 1); refζs=(3, 1, 1, 3)).doit() |> simplify
# 
sum(itr(js)) do two_λs
    amplitude(dc, σs |> StickySymTuple, two_λs; refζs=(3, 1, 1, 3)).doit()^2
end |> simplify
