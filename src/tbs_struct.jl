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

two_J(tbs::ThreeBodySystem) = tbs.two_js[4]
two_j1(tbs::ThreeBodySystem) = tbs.two_js[1]
two_j2(tbs::ThreeBodySystem) = tbs.two_js[2]
two_j3(tbs::ThreeBodySystem) = tbs.two_js[3]

ThreeBodySystem(e,m1,m2,m3; two_jps = ([0,0,0,0], ['+','+','+','+'])) =
    ThreeBodySystem([m1,m2,m3,e], [m1,m2,m3,e].^2,
                    [m2+m3,m3+m1,m1+m2], [m2+m3,m3+m1,m1+m2].^2,
                    [e-m1,e-m2,e-m3], [e-m1,e-m2,e-m3].^2,
                    two_jps[1], two_jps[2]);
#

change_basis_3from1(τ1, tbs::ThreeBodySystem) = change_basis_3from1(τ1..., tbs.msq[1],tbs.msq[2],tbs.msq[3], tbs.msq[4])
change_basis_1from2(τ2, tbs::ThreeBodySystem) = change_basis_3from1(τ2..., tbs.msq[2],tbs.msq[3],tbs.msq[1], tbs.msq[4])
change_basis_2from3(τ3, tbs::ThreeBodySystem) = change_basis_3from1(τ3..., tbs.msq[3],tbs.msq[1],tbs.msq[2], tbs.msq[4])

# coupling scheme
coupling_schemek(k,two_jp,tbs::ThreeBodySystem) = coupling_schemek(k,two_jp,collect(zip(tbs.two_js,tbs.Ps)))
coupling_scheme23(two_jp,tbs::ThreeBodySystem) = coupling_schemek(1,two_jp,tbs)
coupling_scheme12(two_jp,tbs::ThreeBodySystem) = coupling_schemek(3,two_jp,tbs)
coupling_scheme31(two_jp,tbs::ThreeBodySystem) = coupling_schemek(2,two_jp,tbs)

# Dynamic variables
struct DalitzPlotPoint
    σs::SVector{3,Float64}
    two_λs::SVector{4,Int}
end

two_Λ(dpp::DalitzPlotPoint) = dpp.two_λs[4]

DalitzPlotPoint(σ1,σ2,σ3; two_λs::SVector{4,Int}=error("helicities are needed")) = DalitzPlotPoint(SVector(σ1,σ2,σ3), two_λs)
#
DalitzPlotPoint12(σ1,σ2, tbs::ThreeBodySystem; two_λs::SVector{4,Int}=SVector(0,0,0,0)) =
    DalitzPlotPoint(σ1,σ2,gσ3(σ1,σ2,tbs.msq); two_λs = (two_λs == fill(0, 4)) ? tbs.two_js : two_λs) # maximal projection if not provided
DalitzPlotPoint23(σ2,σ3, tbs::ThreeBodySystem; two_λs::SVector{4,Int}=SVector(0,0,0,0)) =
    DalitzPlotPoint(gσ1(σ2,σ3,tbs.msq),σ2,σ3; two_λs = (two_λs == fill(0, 4)) ? tbs.two_js : two_λs) # maximal projection if not provided
DalitzPlotPoint31(σ3,σ1, tbs::ThreeBodySystem; two_λs::SVector{4,Int}=SVector(0,0,0,0)) =
    DalitzPlotPoint(σ1,gσ2(σ3,σ1,tbs.msq),σ3; two_λs = (two_λs == fill(0, 4)) ? tbs.two_js : two_λs) # maximal projection if not provided
#
function randomPoint(tbs::ThreeBodySystem)
    σ1 = tbs.mthsq[1] + rand()*(tbs.sthsq[1]-tbs.mthsq[1])
    σ3 = σ3of1(2rand()-1, σ1, tbs.msq)
    σ2 = gσ2(σ3,σ1,tbs.msq)
    return DalitzPlotPoint(σ1,σ2,σ3; two_λs=SVector([rand(-j:2:j) for j in tbs.two_js]...))
end

function possible_helicities(dpp,tbs)
    [DalitzPlotPoint(dpp.σs...; two_λs = SVector(two_λs...))
        for two_λs = Iterators.product(-tbs.two_js[1]:2:tbs.two_js[1],
                      -tbs.two_js[2]:2:tbs.two_js[2],
                      -tbs.two_js[3]:2:tbs.two_js[3],
                      -tbs.two_js[4]:2:tbs.two_js[4])]
end
