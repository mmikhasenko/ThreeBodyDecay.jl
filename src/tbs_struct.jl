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

ThreeBodySystem(m1,m2,m3,e; two_jps = ([0,0,0,0], ['+','+','+','+'])) =
    e < m1+m2+m3 ? error("Unphysical system, e < sum m_i!") :
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

# dealing with spin 1/2
x2(v) = Int(2v)


#                                                                _|
#    _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|
#  _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|_|_|_|
#  _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|
#    _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|_|
#        _|
#    _|_|

function border(k, tbs; Nx=300)
    (i,j) = ij_from_k(k)
    σiv = range(tbs.mthsq[i], tbs.sthsq[i],length=Nx)
    σkm = [σkofi(k,-1.0,σ,tbs.msq) for σ in σiv]
    σkp = [σkofi(k, 1.0,σ,tbs.msq) for σ in σiv]
    return (σiv, [σkm σkp])
end
#
border31(tbs; Nx=300) = border(3, tbs; Nx=300)
border12(tbs; Nx=300) = border(1, tbs; Nx=300)
border23(tbs; Nx=300) = border(2, tbs; Nx=300)

function flatDalitzPlotSample(k, tbs; Nev::Int=10000, σbins::Int=500)
    (i,j) = ij_from_k(k)
    s = tbs.msq[4]
    density = getbinned1dDensity(σi->sqrt(λ(σi,tbs.msq[j],tbs.msq[k])*λ(σi,s,tbs.msq[i]))/σi, (tbs.mthsq[i],tbs.sthsq[i]), σbins)
    σiv = [rand(density) for _ in 1:Nev]
    σkv = [σkofi(k,2*rand()-1,σ,tbs.msq) for σ in σiv]
    return (σkv, σiv)
end
#
flatDalitzPlotSample31(tbs; Nev::Int=10000, σbins::Int=500) = flatDalitzPlotSample(3, tbs; Nev=Nev, σbins=σbins)
flatDalitzPlotSample12(tbs; Nev::Int=10000, σbins::Int=500) = flatDalitzPlotSample(1, tbs; Nev=Nev, σbins=σbins)
flatDalitzPlotSample23(tbs; Nev::Int=10000, σbins::Int=500) = flatDalitzPlotSample(2, tbs; Nev=Nev, σbins=σbins)
