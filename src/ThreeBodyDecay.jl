module ThreeBodyDecay

using GSL

export gσ1, gσ2, gσ3,
       σ3of1, σ1of2, σ2of3,
       cosθ12, cosθ23, cosθ31,
       cosθhat12,
       cos_mins_θhat2,
       cos_plus_θhat3,
       cosβ,
       λ,
       Kibble, Kibble23, Kibble12, Kibble31

export ThreeBodySystem
export pthsq, pth

# density
export getbinned1dDensity,
       getbinned2dDensity

export flatDalitzPlotSample31

#lineshape
export pole, BW, BWdw

# angual_functions
export ClGd
export jacobi_pols,
       wignerd,
       wignerd_hat
export wignerd_hat_doublearg,
       wignerd_doublearg

include("angle_functions.jl")
include("tbs_struct.jl")
include("rand_corr.jl")
include("rand_gen.jl")
include("lineshape.jl")
include("angular_functions.jl")


end  # module ThreeBodyDecay
