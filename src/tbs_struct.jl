#                                  _|
#    _|_|_|  _|    _|    _|_|_|  _|_|_|_|    _|_|    _|_|_|  _|_|
#  _|_|      _|    _|  _|_|        _|      _|_|_|_|  _|    _|    _|
#      _|_|  _|    _|      _|_|    _|      _|        _|    _|    _|
#  _|_|_|      _|_|_|  _|_|_|        _|_|    _|_|_|  _|    _|    _|
#                  _|
#              _|_|


struct ThreeBodySystem # immutable
    m::SVector{4,Float64}
    msq::SVector{4,Float64}
    #
    mth::SVector{3,Float64}
    mthsq::SVector{3,Float64}
    #
    sth::SVector{3,Float64}
    sthsq::SVector{3,Float64}
    #
    two_js::SVector{4,Int}
    Ps::SVector{4,Char}
end

s(tbs::ThreeBodySystem) = tbs.msq[4]
m0(tbs::ThreeBodySystem) = tbs.m[4]
m1(tbs::ThreeBodySystem) = tbs.m[1]
m2(tbs::ThreeBodySystem) = tbs.m[2]
m3(tbs::ThreeBodySystem) = tbs.m[3]

ThreeBodySystem(e,m1,m2,m3; two_jps = ([0,0,0,0], ['+','+','+','+'])) =
    ThreeBodySystem([m1,m2,m3,e], [m1,m2,m3,e].^2,
                    [m2+m3,m3+m1,m1+m2], [m2+m3,m3+m1,m1+m2].^2,
                    [e-m1,e-m2,e-m3], [e-m1,e-m2,e-m3].^2,
                    two_jps[1], two_jps[2]);
#
# circle(v,i) = (i>0) ? [v[(i+1):end]..., v[1:i]...] : [v[(end+i+1):end]..., v[1:(end+i)]...]
#
# pth(tbs::ThreeBodySystem) = circshift(tbs.m,-1) - circshift(tbs.m,-2)
# pthsq(tbs::ThreeBodySystem) = pth(tbs).^2

change_basis_3from1(τ1, tbs::ThreeBodySystem) = change_basis_3from1(τ1..., tbs.msq[1],tbs.msq[2],tbs.msq[3], tbs.msq[4])
change_basis_1from2(τ2, tbs::ThreeBodySystem) = change_basis_3from1(τ2..., tbs.msq[2],tbs.msq[3],tbs.msq[1], tbs.msq[4])
change_basis_2from3(τ3, tbs::ThreeBodySystem) = change_basis_3from1(τ3..., tbs.msq[3],tbs.msq[1],tbs.msq[2], tbs.msq[4])

#
# Dynamic variables

struct DalitzPlotPoint
    σ123::SVector{3,Float64}
    σ312::SVector{3,Float64}
    σ231::SVector{3,Float64}
end

DalitzPlotPoint(σ1,σ2,σ3) = DalitzPlotPoint(SVector(σ1,σ2,σ3), SVector(σ3,σ1,σ2), SVector(σ2,σ3,σ1))
function randomPoint(tbs::ThreeBodySystem)
    σ1 = tbs.mthsq[1] + rand()*(tbs.sthsq[1]-tbs.mthsq[1])
    σ3 = σ3of1(2rand()-1, σ1, tbs.msq)
    σ2 = gσ2(σ3,σ1,tbs)
    return DalitzPlotPoint(σ1,σ2,σ3)
end

# cosθ23(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s, tbs.msq,                                      dpp.σ123)
# cosθ31(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s,[tbs.msq[2],tbs.msq[3],tbs.msq[1],tbs.msq[4]], dpp.σ231) # (123) permutation
# cosθ12(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθ23(tbs.s,[tbs.msq[3],tbs.msq[1],tbs.msq[2],tbs.msq[4]], dpp.σ312) # (123)^2 permutation

# particle-0 Wigner angle
# cosθhat12(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s,  tbs.msq,                           dpp.σ123)
# cosθhat23(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s, [tbs.msq[2],tbs.msq[3],tbs.msq[1]], dpp.σ231) # (123) permutation
# cosθhat31(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosθhat12(tbs.s, [tbs.msq[3],tbs.msq[1],tbs.msq[2]], dpp.σ312) # (123)^2 permutation

# particle-1 Wigner angle
# cosζ13_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ13_for1(tbs.s,  tbs.msq,                           dpp.σ123)
# cosζ21_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ13_for1(tbs.s, [tbs.msq[1],tbs.msq[3],tbs.msq[2]], SVector(dpp.σ123[1],dpp.σ123[3],dpp.σ123[2]))
# cosζ23_for1(dpp::DalitzPlotPoint,tbs::ThreeBodySystem) = cosζ23_for1(tbs.s,  tbs.msq,                           dpp.σ123)
