#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|


struct ThreeBodySystem
    e::Float64
    s::Float64
    #
    m::Vector{Float64}
    msq::Vector{Float64}
    #
    mth::Vector{Float64}
    mthsq::Vector{Float64}
    #
    sth::Vector{Float64}
    sthsq::Vector{Float64}
end

ThreeBodySystem(e,m1,m2,m3) = ThreeBodySystem(e, e^2, [m1,m2,m3], [m1,m2,m3].^2,
                                            [m2+m3,m3+m1,m1+m2], [m2+m3,m3+m1,m1+m2].^2,
                                            [e-m1,e-m2,e-m3], [e-m1,e-m2,e-m3].^2);
#
circle(v,i) = (i>0) ? [v[(i+1):end]..., v[1:i]...] : [v[(end+i+1):end]..., v[1:(end+i)]...]
#
pth(tbs::ThreeBodySystem) = circle(tbs.m,1) - circle(tbs.m,2)
pthsq(tbs::ThreeBodySystem) = pth(tbs).^2

#
Kibble23(σ2, σ3, tbs::ThreeBodySystem) = Kibble(tbs.s,        tbs.msq,    [tbs.s+sum(tbs.msq)-σ2-σ3,σ2,σ3])
Kibble31(σ3, σ1, tbs::ThreeBodySystem) = Kibble(tbs.s, circle(tbs.msq,1), [σ1,tbs.s+sum(tbs.msq)-σ3-σ1,σ3])
Kibble12(σ1, σ2, tbs::ThreeBodySystem) = Kibble(tbs.s, circle(tbs.msq,2), [σ1,σ2,tbs.s+sum(tbs.msq)-σ1-σ2])
#
σ3of1(σ1,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,       tbs.msq,   σ1,z)
σ1of2(σ3,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,circle(tbs.msq,1),σ2,z)
σ2of3(σ2,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,circle(tbs.msq,2),σ3,z)
#
cosθ12(σ1,σ2,tbs::ThreeBodySystem) = cosθ12(tbs.s,       tbs.msq,          [σ1,σ2,tbs.s+sum(tbs.msq)-σ1-σ2])
cosθ31(σ3,σ1,tbs::ThreeBodySystem) = cosθ12(tbs.s,circle(tbs.msq,1),circle([σ1,tbs.s+sum(tbs.msq)-σ3-σ1,σ3],1))
cosθ23(σ2,σ3,tbs::ThreeBodySystem) = cosθ12(tbs.s,circle(tbs.msq,2),circle([tbs.s+sum(tbs.msq)-σ2-σ3,σ2,σ3],2))
#
