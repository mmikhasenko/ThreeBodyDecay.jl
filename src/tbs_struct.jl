#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|


struct ThreeBodySystem # immutable
    e::Float64
    s::Float64
    #
    m::SVector{3,Float64}
    msq::SVector{3,Float64}
    #
    mth::SVector{3,Float64}
    mthsq::SVector{3,Float64}
    #
    sth::SVector{3,Float64}
    sthsq::SVector{3,Float64}
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
Kibble23(σ2, σ3, tbs::ThreeBodySystem) = Kibble(tbs.s,  tbs.msq,                           [gσ1(σ2,σ3,tbs),σ2,σ3])
Kibble31(σ3, σ1, tbs::ThreeBodySystem) = Kibble(tbs.s, [tbs.msq[2],tbs.msq[3],tbs.msq[1]], [gσ2(σ3,σ1,tbs),σ3,σ1])
Kibble12(σ1, σ2, tbs::ThreeBodySystem) = Kibble(tbs.s, [tbs.msq[3],tbs.msq[1],tbs.msq[2]], [gσ3(σ1,σ2,tbs),σ1,σ2])

#
σ3of1(σ1,z,tbs::ThreeBodySystem) = σ3of1(tbs.s, tbs.msq,                          σ1,z)
σ1of2(σ2,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,[tbs.msq[2],tbs.msq[3],tbs.msq[1]],σ2,z) # (123) permutation
σ2of3(σ3,z,tbs::ThreeBodySystem) = σ3of1(tbs.s,[tbs.msq[3],tbs.msq[1],tbs.msq[2]],σ3,z) # (123)^2 permutation

# Dynamic variables

struct DalitzPlotPoint
    σ123::SVector{3,Float64}
    σ312::SVector{3,Float64}
    σ231::SVector{3,Float64}
end

DalitzPlotPoint(σ1,σ2,σ3) = DalitzPlotPoint(SVector(σ1,σ2,σ3), SVector(σ3,σ1,σ2), SVector(σ2,σ3,σ1))

# scattering angle
cosθ23(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s, tbs.msq,                           dpp.σ123)
cosθ31(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s,[tbs.msq[2],tbs.msq[3],tbs.msq[1]], dpp.σ231) # (123) permutation
cosθ12(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s,[tbs.msq[3],tbs.msq[1],tbs.msq[2]], dpp.σ312) # (123)^2 permutation

# particle-0 Wigner angle
cosθhat12(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s,  tbs.msq,                           dpp.σ123)
cosθhat23(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s, [tbs.msq[2],tbs.msq[3],tbs.msq[1]], dpp.σ231) # (123) permutation
cosθhat31(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s, [tbs.msq[3],tbs.msq[1],tbs.msq[2]], dpp.σ312) # (123)^2 permutation

# particle-1 Wigner angle
cosζ13_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ13_for1(tbs.s,  tbs.msq,                           dpp.σ123)
cosζ21_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ13_for1(tbs.s, [tbs.msq[1],tbs.msq[3],tbs.msq[2]], SVector(dpp.σ123[1],dpp.σ123[3],dpp.σ123[2]))
cosζ23_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ23_for1(tbs.s,  tbs.msq,                           dpp.σ123)
