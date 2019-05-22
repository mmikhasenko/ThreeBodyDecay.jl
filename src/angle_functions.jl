
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

function cosθ23(s,m2s,σs)
    EE4σ = (σs[1]+m2s[2]-m2s[3])*(s-σs[1]-m2s[1])
    pp4σ = sqrt(λ(σs[1],m2s[2],m2s[3])*λ(s,σs[1],m2s[1]))
    rest = σs[3]-m2s[1]-m2s[2]
    return (2σs[1]*rest-EE4σ)/pp4σ
end

function cosθhat12(s,m2s,σs)
    EE4s = (s+m2s[1]-σs[1])*(s+m2s[2]-σs[2])
    pp4s = sqrt(λ(s,m2s[1],σs[1])*λ(s,m2s[2],σs[2]))
    rest = σs[3]-m2s[1]-m2s[2]
    return (EE4s-2s*rest)/pp4s
end

function cosβ(s,m2s,σs)
    EE4m1sq = (s+m2s[1]-σs[1])*(σs[3]-m2s[1]-m2s[2])
    pp4m1sq = sqrt(λ(s,m2s[1],σs[1])*λ(m2s[1],m2s[2],σs[3]))
    rest = σs[2]-s-m2s[2]
    return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
end

function cosθtilde3to1_for1(s,m2s,σs)
    EE4m1sq = (s+m2s[1]-σs[1])*(σs[3]-m2s[1]-m2s[2])
    pp4m1sq = sqrt(λ(s,m2s[1],σs[1])*λ(m2s[1],m2s[2],σs[3]))
    rest = σs[2]-s-m2s[2]
    return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
end

function cosθtilde2to3_for1(s,m2s,σs)
    EE4m1sq = (σs[2]-m2s[3]-m2s[1])*(σs[3]-m2s[1]-m2s[2])
    pp4m1sq = sqrt(λ(σs[2],m2s[3],m2s[1])*λ(σs[3],m2s[1],m2s[2]))
    rest = m2s[2]+m2s[3]-σs[1]
    return (2m2s[1]*rest+EE4m1sq)/pp4m1sq
end
