
#      _|_|                                  _|      _|
#    _|      _|    _|  _|_|_|      _|_|_|  _|_|_|_|        _|_|    _|_|_|      _|_|_|
#  _|_|_|_|  _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|  _|_|
#    _|      _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|      _|_|
#    _|        _|_|_|  _|    _|    _|_|_|      _|_|  _|    _|_|    _|    _|  _|_|_|

λ(x,y,z) = x^2+y^2+z^2 - 2x*y - 2y*z - 2z*x
Kibble(s,m2s,σs) = (2σs[2]*(m2s[1]+s-σs[1])-(s+σs[2]-m2s[2])*(σs[2]+m2s[1]-m2s[3]))^2 - λ(s,σs[2],m2s[2])*λ(σs[2],m2s[1],m2s[3])

function σ3of1(s,m2s,σi,z)
    EE4σ = (σi+m2s[2]-m2s[3])*(s-σi-m2s[1])
    pp4σ = sqrt(λ(σi,m2s[2],m2s[3])*λ(s,σi,m2s[1]))
    return m2s[1]+m2s[2] + (EE4σ + pp4σ*z) / (2σi)
end

# Scattering angle

"""
    Isobar decay angle for the chain-1, i.e. an angle of between \vec p_2 and -\vec p_1 in the (23) rest frame.
"""
function cosθij(k,σs,m2s)
    (i,j) = (k==1 ? (2,3) : (k==2 ? (3,1) : (1,2)))
    #
    s = m2s[4]
    EE4σ = (σs[k]+m2s[i]-m2s[j])*(s-σs[k]-m2s[k])
    pp4σ = sqrt(λ(σs[k],m2s[i],m2s[j])*λ(s,σs[k],m2s[k]))
    rest = σs[j]-m2s[k]-m2s[i]
    return (2σs[k]*rest-EE4σ)/pp4σ
end

cosθ23(s,m2s,σs) = cosθij(1,σs,[m2s...,s])
cosθ31(s,m2s,σs) = cosθij(2,σs,[m2s...,s])
cosθ12(s,m2s,σs) = cosθij(3,σs,[m2s...,s])
# EE4σ = (σs[1]+m2s[2]-m2s[3])*(s-σs[1]-m2s[1])
# pp4σ = sqrt(λ(σs[1],m2s[2],m2s[3])*λ(s,σs[1],m2s[1]))
# rest = σs[3]-m2s[1]-m2s[2]
# return (2σs[1]*rest-EE4σ)/pp4σ

"""
    Wigner angle of the decay particle
"""
# function cosθhat12(s,m2s,σs)
#     EE4s = (s+m2s[1]-σs[1])*(s+m2s[2]-σs[2])
#     pp4s = sqrt(λ(s,m2s[1],σs[1])*λ(s,m2s[2],σs[2]))
#     rest = σs[3]-m2s[1]-m2s[2]
#     return (EE4s-2s*rest)/pp4s
# end
#
function cosθhatk1(k,m2s,σs)
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
cosθhat12(s,m2s,σs) = cosθhatk1(2,[m2s...,s],σs)

# function cosβ(s,m2s,σs)
#     EE4m1sq = (s+m2s[1]-σs[1])*(σs[3]-m2s[1]-m2s[2])
#     pp4m1sq = sqrt(λ(s,m2s[1],σs[1])*λ(m2s[1],m2s[2],σs[3]))
#     rest = σs[2]-s-m2s[2]
#     return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
# end

"""
    Wigner angle of particle 1 to relate chain-1 to chain-3
"""
# function cosζ13_for1(s,m2s,σs)
#     EE4m1sq = (s+m2s[1]-σs[1])*(σs[3]-m2s[1]-m2s[2])
#     pp4m1sq = sqrt(λ(s,m2s[1],σs[1])*λ(m2s[1],m2s[2],σs[3]))
#     rest = σs[2]-s-m2s[2]
#     return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
# end

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
cosζ13_for1(s,m2s,σs) = cosζk1_for1(3,σs,[m2s...,s])
cosζ21_for1(s,m2s,σs) = cosζk1_for1(2,σs,[m2s...,s])

cosζk2_for2(k,σs,m2s) = cosζk1_for1(mod(k-2,3)+1, SVector(σs[2],σs[3],σs[1]), SVector(m2s[2],m2s[3],m2s[1],m2s[4]))
cosζk3_for3(k,σs,m2s) = cosζk1_for1(mod(k,3)+1,   SVector(σs[3],σs[1],σs[2]), SVector(m2s[3],m2s[1],m2s[2],m2s[4]))

"""
    Wigner angle of particle 1 to relate chain-3 to chain-2
"""
function cosζ23_for1(s,m2s,σs)
    EE4m1sq = (σs[2]-m2s[3]-m2s[1])*(σs[3]-m2s[1]-m2s[2])
    pp4m1sq = sqrt(λ(σs[2],m2s[3],m2s[1])*λ(σs[3],m2s[1],m2s[2]))
    rest = m2s[2]+m2s[3]-σs[1]
    return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
end
