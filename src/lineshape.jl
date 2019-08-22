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

function iRho(s,m1sq,m2sq)
    m1, m2 = sqrt(m1sq), sqrt(m2sq)
    #
    sth,spth = (m1+m2)^2, (m1-m2)^2;
    λh = sqrt(s-sth)*sqrt(s-spth)
    return 1im*λh/s
end

#  _|  _|                                _|
#  _|      _|_|_|      _|_|      _|_|_|  _|_|_|      _|_|_|  _|_|_|      _|_|
#  _|  _|  _|    _|  _|_|_|_|  _|_|      _|    _|  _|    _|  _|    _|  _|_|_|_|
#  _|  _|  _|    _|  _|            _|_|  _|    _|  _|    _|  _|    _|  _|
#  _|  _|  _|    _|    _|_|_|  _|_|_|    _|    _|    _|_|_|  _|_|_|      _|_|_|
#                                                            _|
#                                                            _|

pole(σ,mcsq) = 1.0/(mcsq - σ)
BW(σ,m,Γ) = m*Γ*pole(σ,m^2-1im*m*Γ)

function BWdw(σ,m,Γ,m1,m2)
    √σ < (m1+m2) && return 0.0im
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
