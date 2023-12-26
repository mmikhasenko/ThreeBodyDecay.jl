#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|

const MassTuple{T} = NamedTuple{(:m1, :m2, :m3, :m0),NTuple{4,T}}
function ThreeBodyMasses(m1, m2, m3; m0)
    tm0 = typeof(m0)
    tm0 <: Number && (m0 < m1 + m2 + m3) && error("m₀ should be bigger than m₁+m₂+m₃")
    MassTuple{tm0}((m1, m2, m3, m0))
end
# 
ThreeBodyMasses(; m1, m2, m3, m0) = ThreeBodyMasses(m1, m2, m3; m0)

lims(k::Int, ms::MassTuple) = k == 1 ? lims1(ms) : ((k == 2) ? lims2(ms) : lims3(ms))
lims1(ms::MassTuple) = ((ms.m2 + ms.m3)^2, (ms.m0 - ms.m1)^2)
lims2(ms::MassTuple) = ((ms.m3 + ms.m1)^2, (ms.m0 - ms.m2)^2)
lims3(ms::MassTuple) = ((ms.m1 + ms.m2)^2, (ms.m0 - ms.m3)^2)
#
import Base: getindex, ^, length, iterate
^(ms::MassTuple, i::Int) = Tuple(ms) .^ i

# -----------------------------------------------------
const SpinTuple = NamedTuple{(:two_h1, :two_h2, :two_h3, :two_h0),NTuple{4,Int}}

ThreeBodySpins(two_h1, two_h2, two_h3; two_h0=error("used the format ThreeBodySpins(1,1,0; two_h0=2)")) =
    isodd(two_h1 + two_h2 + two_h3 + two_h0) ?
    error("baryon number is not conserved") :
    SpinTuple((two_h1, two_h2, two_h3, two_h0))
# 
# 
@with_kw struct ThreeBodySystem{T,K}
    ms::T
    two_js::K = ThreeBodySpins(0, 0, 0; two_h0=0)
end
#
# convenient constructors
ThreeBodySystem(ms::MassTuple) = ThreeBodySystem(ms=ms)
ThreeBodySystem(m1, m2, m3; m0, two_js=ThreeBodySpins(0, 0, 0; two_h0=0)) = ThreeBodySystem(ThreeBodyMasses(m1, m2, m3; m0=m0), two_js)
#
two_j0(tbs::ThreeBodySystem) = tbs.two_js[4]
two_j1(tbs::ThreeBodySystem) = tbs.two_js[1]
two_j2(tbs::ThreeBodySystem) = tbs.two_js[2]
two_j3(tbs::ThreeBodySystem) = tbs.two_js[3]

# -----------------------------------------------------

const ParityTuple = NamedTuple{(:P1, :P2, :P3, :P0),NTuple{4,Char}}
#
ThreeBodyParities(P1, P2, P3;
    P0=error("used the format ThreeBodyParities('+','-','+'; P0='±')")) =
    ParityTuple((P1, P2, P3, P0))
# -----------------------------------------------------

# Dynamic variables
const MandestamTuple{T} = NamedTuple{(:σ1, :σ2, :σ3),NTuple{3,T}}

"""
    Invariants(ms::MassTuple{T}; σ1, σ2)
    Invariants(ms::MassTuple{T}; σ1, σ3)
    Invariants(ms::MassTuple{T}; σ2, σ3)

Construct a tuple of (σ1, σ2, σ3) from just two invariants and the mass tuple.
"""
function Invariants(ms::MassTuple{T};
    σ1=-one(ms.m0), σ2=-one(ms.m0), σ3=-one(ms.m0)) where {T}
    # 
    !((σ1 == -one(ms.m0)) || (σ2 == -one(ms.m0)) || (σ3 == -one(ms.m0))) &&
        error("the method works with TWO invariants given: $((σ1,σ2,σ3))")
    # 
    σ3 == -one(ms.m0) && return MandestamTuple{T}((σ1, σ2, σ3=sum(ms^2) - σ1 - σ2))
    σ1 == -one(ms.m0) && return MandestamTuple{T}((sum(ms^2) - σ2 - σ3, σ2, σ3))
    return MandestamTuple{T}((σ1=σ1, σ2=sum(ms^2) - σ3 - σ1, σ3=σ3))
end
Invariants(; σ1, σ2, σ3) = MandestamTuple{typeof(σ1)}((σ1, σ2, σ3))
Invariants(σ1, σ2, σ3) = MandestamTuple{typeof(σ1)}((σ1, σ2, σ3))

# -----------------------------------------------------

@with_kw struct DalitzPlotPoint{I,S}
    σs::I
    two_λs::S
end
#
function randomPoint(ms::MassTuple)
    σ1 = lims1(ms)[1] + rand() * (lims1(ms)[2] - lims1(ms)[1])
    σ3 = σ3of1(2rand() - 1, σ1, ms^2)
    return Invariants(ms; σ1=σ1, σ3=σ3)
end

function randomPoint(tbs::ThreeBodySystem)
    DalitzPlotPoint(σs=randomPoint(tbs.ms),
        two_λs=[rand(-j:2:j) for j in tbs.two_js])
end

#                                                                _|
#    _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|
#  _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|_|_|_|
#  _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|
#    _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|_|
#        _|
#    _|_|

"""
    polardalitz2invariants(θ, expansion_point)

For given polar angle θ, it returns an (σ1,σ2,σ3) Tuple of polynomials of radius r(θ) around the expansion point.
The polynomial works as a function of the r coordinate.
"""
polardalitz2invariants(θ, expansion_point::Tuple) =
    (Polynomial([0, -cos(θ)]),
        Polynomial([0, cos(θ + π / 3)]),
        Polynomial([0, cos(θ - π / 3)])) .+ expansion_point

function border(ms::MassTuple{T}; Nx::Int=300) where {T}
    # 
    expansion_point = let
        f = 0.5
        z = 0.0
        σ1 = (ms[2] + ms[3])^2 + f * ((ms[4] - ms[1])^2 - (ms[2] + ms[3])^2)
        σ3 = σ3of1(z, σ1, ms^2)
        Invariants(ms; σ1, σ3)
    end
    # 
    σs(θ) = polardalitz2invariants(θ, expansion_point |> Tuple)
    ϕ0 = Kibble(expansion_point, ms^2)
    ϕ(σs) = Kibble(σs, ms^2)
    # 
    function rborder(θ)
        _roots = PolynomialRoots.roots(coeffs(ϕ(σs(θ))))
        filter(_roots) do r
            (abs(imag(r)) < 1e-10) && real(r) > 0.0
        end |> real |> minimum
    end
    function σsborder(θ) # evaluate the polynomials
        r = rborder(θ)
        return map(P -> P(r), σs(θ))
    end
    θs = range(-π / 9, 2π - π / 9, length=Nx)
    σs_tuple = σsborder.(θs)
    return MandestamTuple{T}.(σs_tuple)
end

# border13, border12, border21, border23, border32
for (i, j) in ((1, 2), (2, 1), (2, 3), (3, 2), (3, 1), (1, 3))
    eval(quote
        $(Symbol(:border, i, j))(ms; Nx::Int=300) =
            NamedTuple{$(Symbol(:σ, i), Symbol(:σ, j))}.(border(ms; Nx))
    end)
end

# 
function flatDalitzPlotSample(ms::MassTuple; Nev::Int=10000, σbins::Int=500)
    @unpack m0, m1, m2, m3 = ms
    density = getbinned1dDensity(σ1 -> sqrt(Kallen(σ1, m2^2, m3^2) * Kallen(σ1, m0^2, m1^2)) / σ1, lims1(ms), σbins)
    σ1v = [rand(density) for _ in 1:Nev]
    σ3v = [σ3of1(2 * rand() - 1, σ1, ms^2) for σ1 in σ1v]
    return [Invariants(ms; σ1=σ1, σ3=σ3) for (σ1, σ3) in zip(σ1v, σ3v)]
end

#
inrange(x, r) = r[1] < x < r[2]
inphrange(σs::MandestamTuple, ms::MassTuple) = Kibble(σs, ms^2) < 0 &&
                                               inrange(σs[1], lims1(ms)) && inrange(σs[2], lims2(ms)) && inrange(σs[3], lims3(ms))
#
# 
change_basis_3from1(τ1, ms::MassTuple) = change_basis_3from1(τ1..., ms.m1^2, ms.m2^2, ms.m3^2, ms.m0^2)
change_basis_1from2(τ2, ms::MassTuple) = change_basis_3from1(τ2..., ms.m2^2, ms.m3^2, ms.m1^2, ms.m0^2)
change_basis_2from3(τ3, ms::MassTuple) = change_basis_3from1(τ3..., ms.m3^2, ms.m1^2, ms.m2^2, ms.m0^2)


