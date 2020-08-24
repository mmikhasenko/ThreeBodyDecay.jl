
#      _|_|                                  _|      _|
#    _|      _|    _|  _|_|_|      _|_|_|  _|_|_|_|        _|_|    _|_|_|      _|_|_|
#  _|_|_|_|  _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|  _|_|
#    _|      _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|      _|_|
#    _|        _|_|_|  _|    _|    _|_|_|      _|_|  _|    _|_|    _|    _|  _|_|_|

λ(x,y,z) = x^2+y^2+z^2 - 2x*y - 2y*z - 2z*x
Kibble(σs,msq) = λ(λ(msq[4],msq[1],σs[1]),
                   λ(msq[4],msq[2],σs[2]),
                   λ(msq[4],msq[3],σs[3]))
#
ij_from_k(k) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
#
function σkofi(k,z,σi,msq)
    (i,j) = ij_from_k(k)
    #
    s = msq[4]
    EE4σ = (σi+msq[j]-msq[k])*(s-σi-msq[i])
    pp4σ_sq = λ(σi,msq[j],msq[k])*λ(s,σi,msq[i])
    pp4σ_sq = (pp4σ_sq < 0) ? 0.0 : pp4σ_sq # for numerical errors
    σk = msq[i]+msq[j] + (EE4σ + sqrt(pp4σ_sq)*z) / (2σi)
    return σk
end
σ3of1(z,σ1,msq) = σkofi(3,z,σ1,msq)
σ1of2(z,σ2,msq) = σkofi(1,z,σ2,msq)
σ2of3(z,σ3,msq) = σkofi(2,z,σ3,msq)

# Scattering angle

"""
    Isobar decay angle for the chain-1, i.e. an angle of between \vec p_2 and -\vec p_1 in the (23) rest frame.
"""
function cosθij(k,σs,msq)
    (i,j) = ij_from_k(k)
    #
    s = msq[4]
    EE4σ = (σs[k]+msq[i]-msq[j])*(s-σs[k]-msq[k])
    pp4σ = sqrt(λ(σs[k],msq[i],msq[j])*λ(s,σs[k],msq[k]))
    rest = σs[j]-msq[k]-msq[i]
    return (2σs[k]*rest-EE4σ)/pp4σ
end
cosθ23(σs,msq) = cosθij(1,σs,msq)
cosθ31(σs,msq) = cosθij(2,σs,msq)
cosθ12(σs,msq) = cosθij(3,σs,msq)

"""
    Wigner angle of the decay particle
"""
function cosθhatk1(k,σs,msq)
    k==1 && return 1.0
    #
    (i,j) = (k==2 ? (1,3) : (1,2)) # k=2 are flipped coz k1=21 is negative
    #
    s = msq[4]
    EE4s = (s+msq[i]-σs[i])*(s+msq[k]-σs[k])
    pp4s = sqrt(λ(s,msq[i],σs[i])*λ(s,msq[k],σs[k]))
    rest = σs[j]-msq[i]-msq[k]
    return (EE4s-2s*rest)/pp4s
end
cosθhatk2(k,σs,msq) = cosθhatk1(mod(k-2,3)+1, SVector(σs[2],σs[3],σs[1]), SVector(msq[2],msq[3],msq[1],msq[4]))
cosθhatk3(k,σs,msq) = cosθhatk1(mod(k,3)+1,   SVector(σs[3],σs[1],σs[2]), SVector(msq[3],msq[1],msq[2],msq[4]))
#
cosθhat31(σs,msq) = cosθhatk1(3,σs,msq)
cosθhat12(σs,msq) = cosθhatk1(2,σs,msq)
cosθhat23(σs,msq) = cosθhatk2(3,σs,msq)

"""
    Wigner angle of particle 1 to relate chain-1 to chain-3
"""

function cosζk1_for1(k,σs,msq)
    k==1 && return 1.0
    #
    (i,j) = (k==2 ? (1,3) : (1,2)) # k=2 are flipped coz k1=21 is negative
    #
    s = msq[4]
    EE4m1sq = (s+msq[i]-σs[i])*(σs[k]-msq[i]-msq[j])
    pp4m1sq = sqrt(λ(s,msq[i],σs[i])*λ(msq[i],msq[j],σs[k]))
    rest = σs[j]-s-msq[j]
    return (2msq[i]*rest+EE4m1sq)/pp4m1sq
end
cosζk2_for2(k,σs,msq) = cosζk1_for1(mod(k-2,3)+1, SVector(σs[2],σs[3],σs[1]), SVector(msq[2],msq[3],msq[1],msq[4]))
cosζk3_for3(k,σs,msq) = cosζk1_for1(mod(k,3)+1,   SVector(σs[3],σs[1],σs[2]), SVector(msq[3],msq[1],msq[2],msq[4]))
# particular
cosζ13_for1(σs,msq) = cosζk1_for1(3,σs,msq)
cosζ21_for1(σs,msq) = cosζk1_for1(2,σs,msq)
cosζ21_for2(σs,msq) = cosζk2_for2(1,σs,msq)
cosζ32_for2(σs,msq) = cosζk2_for2(3,σs,msq)
cosζ32_for3(σs,msq) = cosζk2_for2(2,σs,msq)
cosζ13_for3(σs,msq) = cosζk2_for2(1,σs,msq)

"""
    Wigner angle of particle 1 to relate chain-3 to chain-2
"""
function cosζ23_for1(σs,msq)
    msq[1] ≈ 0 && return 1.0
    s = msq[4]
    EE4m1sq = (σs[2]-msq[3]-msq[1])*(σs[3]-msq[1]-msq[2])
    pp4m1sq = sqrt(λ(σs[2],msq[3],msq[1])*λ(σs[3],msq[1],msq[2]))
    rest = msq[2]+msq[3]-σs[1]
    return (2msq[1]*rest+EE4m1sq)/pp4m1sq
end
cosζ31_for2(σs,msq) = cosζ23_for1(SVector(σs[2],σs[3],σs[1]), SVector(msq[2],msq[3],msq[1],msq[4]))
cosζ12_for3(σs,msq) = cosζ23_for1(SVector(σs[3],σs[1],σs[2]), SVector(msq[3],msq[1],msq[2],msq[4]))

"""
    Phase for wigner d-functions for clockwise rotations
"""
phase(two_λ1_minus_λ2) = (abs(two_λ1_minus_λ2) % 4 == 2 ? -1.0 : 1.0)
phase(two_λ1,two_λ2) = phase(two_λ1 - two_λ2)

"""
    Phase for wigner d-functions for clockwise rotations
        with a check if indices are in the sequential order.
"""
function phase(i,j,two_λ1,two_λ2)
    (i==j) && return (two_λ1 == two_λ2 ? 1.0 : 0.0)
    ((i,j)==(1,2) || (i,j)==(2,3) || (i,j)==(3,1)) && return 1.0
    return phase(two_λ1,two_λ2)
end

"""
    Calculate normalized values in square coordinates
"""
function squaredalitz(k,σs,tbs)
    cθ = cosθij(k,σs,tbs.msq)
    cθ = cθ ≥ 1.0 ? 1.0 : (cθ ≤ -1.0 ? -1.0 : cθ)
    y = acos(cθ) / π
    xn = 2*(sqrt(σs[k]) - tbs.mth[k]) / (tbs.sth[k]-tbs.mth[k])-1
    x = acos(xn) / π
    return (x=x, y=y)
end

squaredalitz1(σs,tbs) = squaredalitz(1,σs,tbs)
squaredalitz2(σs,tbs) = squaredalitz(2,σs,tbs)
squaredalitz3(σs,tbs) = squaredalitz(3,σs,tbs)
