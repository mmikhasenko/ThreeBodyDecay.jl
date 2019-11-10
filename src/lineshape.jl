function ChewMandestam(s,m1sq,m2sq)
    m1, m2 = sqrt(m1sq), sqrt(m2sq)
    #
    sth,spth = (m1+m2)^2, (m1-m2)^2;
    λh = sqrt(spth-s)*sqrt(sth-s)
    #
    val = 1/(π) * (
        λh/s*log((m1sq+m2sq-s+λh)/(2*m1*m2))-
            (m1sq-m2sq) * (1.0/s - 1/sth) * log(m1/m2)
    )
    return val
end

function Rho(s,m1sq,m2sq)
    m1, m2 = sqrt(m1sq), sqrt(m2sq)
    #
    sth,spth = (m1+m2)^2, (m1-m2)^2;
    λh = sqrt(s-sth)*sqrt(s-spth)
    return λh/s
end
iRho(s,m1sq,m2sq) = 1im*Rho(s,m1sq,m2sq)


#  _|  _|                                _|
#  _|      _|_|_|      _|_|      _|_|_|  _|_|_|      _|_|_|  _|_|_|      _|_|
#  _|  _|  _|    _|  _|_|_|_|  _|_|      _|    _|  _|    _|  _|    _|  _|_|_|_|
#  _|  _|  _|    _|  _|            _|_|  _|    _|  _|    _|  _|    _|  _|
#  _|  _|  _|    _|    _|_|_|  _|_|_|    _|    _|    _|_|_|  _|_|_|      _|_|_|
#                                                            _|
#                                                            _|

pole(σ,mcsq) = 1.0/(mcsq - σ)
BW(σ,m,Γ) = pole(σ,m^2-1im*m*Γ)

function BWdw(σ,m,Γ,m1,m2)
    √σ <= (m1+m2) && return 0.0im
    ρσ = sqrt(λ(σ,m1^2,m2^2))/σ
    ρ0 = sqrt(λ(m^2,m1^2,m2^2))/m^2
    pole(σ,m^2-1im*m*Γ*ρσ/ρ0)
end

struct Lineshape{T<:Real}
    itype::Int
    type::String
    pars::Vector{T}
end

# Breit-Wigner with constant width
Lineshape(m,Γ) = Lineshape(0,"BW",[m,Γ])
BreitWigner(m,Γ) = Lineshape(0,"BW",[m,Γ])

function amp(s,lsh::Lineshape)
    (lsh.itype == 0) && return BW(s,lsh.pars[1],lsh.pars[2])
    error("No itype #$(lsh.itype) found")
    return 0.0im
end


#            _|
#  _|_|_|    _|_|_|          _|_|_|  _|_|_|
#  _|    _|  _|    _|      _|_|      _|    _|
#  _|    _|  _|    _|          _|_|  _|    _|
#  _|_|_|    _|    _|  _|  _|_|_|    _|_|_|    _|
#  _|                                _|
#  _|                                _|

function RhoQTB(s,Xlineshape,msq; channel::Int=1)
    m1sq, m2sq, m3sq = msq
    m1, m2, m3 = sqrt.(msq)
    #
    k = channel
    i,j = ij_from_k(k)
    #
    (√s < (m1+m2+m3) || √s ≈ (m1+m2+m3)) && return 0.0
    val = quadgk(σ->abs2(Xlineshape(s,σ))*sqrt(λ(s,σ,msq[k])*λ(σ,msq[i],msq[j]))/σ, (√msq[i]+√msq[j])^2, (√s-√msq[k])^2)[1]
    return val / s # /(2π*(8π)^2)
end

RhoQTB(s,m1,m2,Γ1,m1th) = RhoQTB(s,(s,σ)->BW(σ,m1,Γ1),[m2,m1th/2,m1th/2].^2)
iRhoQTB(s,m1,m2,Γ1,m1th) = 1im*RhoQTB(s,m1,m2,Γ1,m1th)

function three_body_phase_space_integral(function_σs, tbs)
    msq = tbs.msq
    m1sq,m2sq,m3sq,s = msq
    #
    σ1min,σ1max = tbs.mthsq[1],tbs.sthsq[1]
    function integrand(x,f)
        σ1 = σ1min + x[1]*(σ1max-σ1min)
        z = 2*x[2]-1
        σ3 = σ3of1(z,σ1,msq); σ2 = gσ2(σ3,σ1,msq)
        σs = [σ1, σ2, σ3]
        #
        r = function_σs(σs) *
            sqrt(λ(s,σ1, m1sq)*λ(σ1,m2sq,m3sq))/σ1 *
            (σ1max-σ1min) # jacobians
        f[1],f[2] = reim(r)
    end
    return complex(cuhre(integrand,2,2)[1]...) / s # /(2π*(8π)^2)
end
