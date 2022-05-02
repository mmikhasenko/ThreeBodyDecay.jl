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
    ρσ = sqrt(Kallen(σ,m1^2,m2^2))/σ
    ρ0 = sqrt(Kallen(m^2,m1^2,m2^2))/m^2
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
    (lsh.itype == 9) && return ScattLenT11(s,lsh.pars[1:3],lsh.pars[4:5],lsh.pars[6:7])
    error("No itype #$(lsh.itype) found")
    return 0.0im
end

# Scattering length approximation
ScattLenApproximation(;
    Ms::Tuple{Real,Real,Real} = error("give Ms = (c11, c12, c22)"),
    ms1::Tuple{Real,Real} = error("give ms1 = (m1, m2)"),
    ms2::Tuple{Real,Real}) = Lineshape(9,"ScattLen",[Ms...,ms1...,ms2...])

function ScattLenMatrix(s,Ms,ms1,ms2)
    M  = SMatrix{2,2}(Ms[1],Ms[2],Ms[2],Ms[3])
    iρ = SMatrix{2,2}(-sqrt((ms1[1]+ms1[2])^2-s)*sqrt(s-(ms1[1]-ms1[2])^2)/s, 0.0im, 0.0im,
                      -sqrt((ms2[1]+ms2[2])^2-s)*sqrt(s-(ms2[1]-ms2[2])^2)/s)
    return inv(M - iρ)
end
ScattLenT11(s,Ms,ms1,ms2) = ScattLenMatrix(s,Ms,ms1,ms2)[1,1]

#            _|
#  _|_|_|    _|_|_|          _|_|_|  _|_|_|
#  _|    _|  _|    _|      _|_|      _|    _|
#  _|    _|  _|    _|          _|_|  _|    _|
#  _|_|_|    _|    _|  _|  _|_|_|    _|_|_|    _|
#  _|                                _|
#  _|                                _|

function RhoQTB(s,Xlineshape,msq; channel::Int=1)
    m1, m2, m3 = sqrt.(msq)
    #
    k = channel
    i,j = ij_from_k(k)
    #
    (√s < (m1+m2+m3) || √s ≈ (m1+m2+m3)) && return 0.0
    val = quadgk(σ->abs2(Xlineshape(σ))*sqrt(Kallen(s,σ,msq[k])*Kallen(σ,msq[i],msq[j]))/σ, (√msq[i]+√msq[j])^2, (√s-√msq[k])^2)[1]
    return val / s # /(2π*(8π)^2)
end

RhoQTB(s,m1,m2,Γ1,m1th) = RhoQTB(s,σ->BW(σ,m1,Γ1),[m2,m1th/2,m1th/2].^2)
iRhoQTB(s,m1,m2,Γ1,m1th) = 1im*RhoQTB(s,m1,m2,Γ1,m1th)

function three_body_phase_space_integral(function_σs, ms)
    m1sq,m2sq,m3sq,s = ms^2
    #
    σ1min,σ1max = lims1(ms)
    function integrand(x,f)
        σ1 = σ1min + x[1]*(σ1max-σ1min)
        z = 2*x[2]-1
        σs = Invariants(ms;σ1=σ1,σ3=σ3of1(z,σ1,ms^2))
        #
        r = function_σs(σs) *
            sqrt(Kallen(s,σ1, m1sq)*Kallen(σ1,m2sq,m3sq))/σ1 *
            (σ1max-σ1min) # jacobians
        f[1],f[2] = reim(r)
    end
    return complex(cuhre(integrand,2,2)[1]...) / s # /(2π*(8π)^2)
end
