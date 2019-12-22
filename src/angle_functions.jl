
#      _|_|                                  _|      _|
#    _|      _|    _|  _|_|_|      _|_|_|  _|_|_|_|        _|_|    _|_|_|      _|_|_|
#  _|_|_|_|  _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|  _|_|
#    _|      _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|      _|_|
#    _|        _|_|_|  _|    _|    _|_|_|      _|_|  _|    _|_|    _|    _|  _|_|_|

gσ1(σ2,σ3,msq) = sum(msq)-σ2-σ3
gσ2(σ3,σ1,msq) = sum(msq)-σ3-σ1
gσ3(σ1,σ2,msq) = sum(msq)-σ1-σ2

#

λ(x,y,z) = x^2+y^2+z^2 - 2x*y - 2y*z - 2z*x
Kibble(σs,m2s) = (2σs[2]*(m2s[1]+m2s[4]-σs[1])-(m2s[4]+σs[2]-m2s[2])*(σs[2]+m2s[1]-m2s[3]))^2 - λ(m2s[4],σs[2],m2s[2])*λ(σs[2],m2s[1],m2s[3])
#
Kibble23(σ2, σ3, msq) = Kibble(SVector(gσ1(σ2,σ3,msq),σ2,σ3), msq)
Kibble31(σ3, σ1, msq) = Kibble(SVector(gσ2(σ3,σ1,msq),σ3,σ1), SVector(msq[2],msq[3],msq[1],msq[4]))
Kibble12(σ1, σ2, msq) = Kibble(SVector(gσ3(σ1,σ2,msq),σ1,σ2), SVector(msq[3],msq[1],msq[2],msq[4]))
#
ij_from_k(k) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
#
function σkofi(k,z,σi,m2s)
    (i,j) = ij_from_k(k)
    #
    s = m2s[4]
    EE4σ = (σi+m2s[j]-m2s[k])*(s-σi-m2s[i])
    pp4σ_sq = λ(σi,m2s[j],m2s[k])*λ(s,σi,m2s[i])
    pp4σ_sq = (pp4σ_sq < 0) ? 0.0 : pp4σ_sq # for numerical errors
    σk = m2s[i]+m2s[j] + (EE4σ + sqrt(pp4σ_sq)*z) / (2σi)
    return σk
end
σ3of1(z,σ1,m2s) = σkofi(3,z,σ1,m2s)
σ1of2(z,σ2,m2s) = σkofi(1,z,σ2,m2s)
σ2of3(z,σ3,m2s) = σkofi(2,z,σ3,m2s)

# Scattering angle

"""
    Isobar decay angle for the chain-1, i.e. an angle of between \vec p_2 and -\vec p_1 in the (23) rest frame.
"""
function cosθij(k,σs,m2s)
    (i,j) = ij_from_k(k)
    #
    s = m2s[4]
    EE4σ = (σs[k]+m2s[i]-m2s[j])*(s-σs[k]-m2s[k])
    pp4σ = sqrt(λ(σs[k],m2s[i],m2s[j])*λ(s,σs[k],m2s[k]))
    rest = σs[j]-m2s[k]-m2s[i]
    return (2σs[k]*rest-EE4σ)/pp4σ
end
cosθ23(σs,m2s) = cosθij(1,σs,m2s)
cosθ31(σs,m2s) = cosθij(2,σs,m2s)
cosθ12(σs,m2s) = cosθij(3,σs,m2s)

"""
    Wigner angle of the decay particle
"""
function cosθhatk1(k,σs,m2s)
    k==1 && return 1.0
    #
    (i,j) = (k==2 ? (1,3) : (1,2)) # k=2 are flipped coz k1=21 is negative
    #
    s = m2s[4]
    EE4s = (s+m2s[i]-σs[i])*(s+m2s[k]-σs[k])
    pp4s = sqrt(λ(s,m2s[i],σs[i])*λ(s,m2s[k],σs[k]))
    rest = σs[j]-m2s[i]-m2s[k]
    return (EE4s-2s*rest)/pp4s
end
cosθhatk2(k,σs,m2s) = cosθhatk1(mod(k-2,3)+1, SVector(σs[2],σs[3],σs[1]), SVector(m2s[2],m2s[3],m2s[1],m2s[4]))
cosθhatk3(k,σs,m2s) = cosθhatk1(mod(k,3)+1,   SVector(σs[3],σs[1],σs[2]), SVector(m2s[3],m2s[1],m2s[2],m2s[4]))
#
cosθhat31(σs,m2s) = cosθhatk1(3,σs,m2s)
cosθhat12(σs,m2s) = cosθhatk1(2,σs,m2s)
cosθhat23(σs,m2s) = cosθhatk2(3,σs,m2s)

"""
    Wigner angle of particle 1 to relate chain-1 to chain-3
"""

function cosζk1_for1(k,σs,m2s)
    k==1 && return 1.0
    #
    (i,j) = (k==2 ? (1,3) : (1,2)) # k=2 are flipped coz k1=21 is negative
    #
    s = m2s[4]
    EE4m1sq = (s+m2s[i]-σs[i])*(σs[k]-m2s[i]-m2s[j])
    pp4m1sq = sqrt(λ(s,m2s[i],σs[i])*λ(m2s[i],m2s[j],σs[k]))
    rest = σs[j]-s-m2s[j]
    return (2m2s[i]*rest+EE4m1sq)/pp4m1sq
end
cosζk2_for2(k,σs,m2s) = cosζk1_for1(mod(k-2,3)+1, SVector(σs[2],σs[3],σs[1]), SVector(m2s[2],m2s[3],m2s[1],m2s[4]))
cosζk3_for3(k,σs,m2s) = cosζk1_for1(mod(k,3)+1,   SVector(σs[3],σs[1],σs[2]), SVector(m2s[3],m2s[1],m2s[2],m2s[4]))
# particular
cosζ13_for1(σs,m2s) = cosζk1_for1(3,σs,m2s)
cosζ21_for1(σs,m2s) = cosζk1_for1(2,σs,m2s)
cosζ21_for2(σs,m2s) = cosζk2_for2(1,σs,m2s)
cosζ32_for2(σs,m2s) = cosζk2_for2(3,σs,m2s)
cosζ32_for3(σs,m2s) = cosζk2_for2(2,σs,m2s)
cosζ13_for3(σs,m2s) = cosζk2_for2(1,σs,m2s)

"""
    Wigner angle of particle 1 to relate chain-3 to chain-2
"""
function cosζ23_for1(σs,m2s)
    m2s[1] ≈ 0 && return 1.0
    s = m2s[4]
    EE4m1sq = (σs[2]-m2s[3]-m2s[1])*(σs[3]-m2s[1]-m2s[2])
    pp4m1sq = sqrt(λ(σs[2],m2s[3],m2s[1])*λ(σs[3],m2s[1],m2s[2]))
    rest = m2s[2]+m2s[3]-σs[1]
    return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
end
cosζ31_for2(σs,m2s) = cosζ23_for1(SVector(σs[2],σs[3],σs[1]), SVector(m2s[2],m2s[3],m2s[1],m2s[4]))
cosζ12_for3(σs,m2s) = cosζ23_for1(SVector(σs[3],σs[1],σs[2]), SVector(m2s[3],m2s[1],m2s[2],m2s[4]))

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
