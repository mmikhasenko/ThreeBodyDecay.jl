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
# circle(v,i) = (i>0) ? [v[(i+1):end]..., v[1:i]...] : [v[(end+i+1):end]..., v[1:(end+i)]...]
#
pth(tbs::ThreeBodySystem) = circshift(tbs.m,-1) - circshift(tbs.m,-2)
pthsq(tbs::ThreeBodySystem) = pth(tbs).^2

gσ1(σ2,σ3,tbs) = tbs.s+sum(tbs.msq)-σ2-σ3
gσ2(σ3,σ1,tbs) = tbs.s+sum(tbs.msq)-σ3-σ1
gσ3(σ1,σ2,tbs) = tbs.s+sum(tbs.msq)-σ1-σ2
#
Kibble23(σ2, σ3, tbs::ThreeBodySystem) = Kibble(tbs.s,           tbs.msq,              [σ1,σ2,gσ3(σ1,σ2,tbs)]   )
Kibble31(σ3, σ1, tbs::ThreeBodySystem) = Kibble(tbs.s, circshift(tbs.msq,-1), circshift([σ1,gσ2(σ3,σ1,tbs),σ3],-1))
Kibble12(σ1, σ2, tbs::ThreeBodySystem) = Kibble(tbs.s, circshift(tbs.msq,-2), circshift([gσ1(σ2,σ3,tbs),σ2,σ3],-2))
#
σ3of1(σ1,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,          tbs.msq,   σ1,z)
σ1of2(σ3,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,circshift(tbs.msq,-1),σ3,z)
σ2of3(σ2,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,circshift(tbs.msq,-2),σ2,z)
#
cosθ12(σ1,σ2,tbs::ThreeBodySystem) = cosθ12(tbs.s,          tbs.msq,              [σ1,σ2,gσ3(σ1,σ2,tbs)]   )
cosθ31(σ3,σ1,tbs::ThreeBodySystem) = cosθ12(tbs.s,circshift(tbs.msq,-1), circshift([σ1,gσ2(σ3,σ1,tbs),σ3],-1))
cosθ23(σ2,σ3,tbs::ThreeBodySystem) = cosθ12(tbs.s,circshift(tbs.msq,-2), circshift([gσ1(σ2,σ3,tbs),σ2,σ3],-2))
#
cos_mins_θhat2(σ3,σ1,tbs::ThreeBodySystem) = cosθhat12(tbs.s,           tbs.msq    ,           [σ1,gσ2(σ3,σ1,tbs),σ3])
cos_plus_θhat3(σ1,σ2,tbs::ThreeBodySystem) = cosθhat12(tbs.s, circshift(tbs.msq,1), circshift([σ1,σ2,gσ3(σ1,σ2,tbs)],1))
