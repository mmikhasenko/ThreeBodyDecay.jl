module ThreeBodyDecay

using GSL

export gσ1, gσ2, gσ3,
       σ3of1, σ1of2, σ2of3,
       cosθ12, cosθ23, cosθ31,
       cosθhat12,
       cosθhat23,
       cosθhat31,
       cos_plus_θhat31,
       cos_mins_θhat21,
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
       wignerD,
       wignerd_hat
export wignerd_hat_doublearg,
       wignerd_doublearg,
       wignerD_doublearg

 export change_basis_3from1,
        change_basis_1from2,
        change_basis_2from3

export cosθtilde3to1_for1,
       cosθtilde2to3_for1,
       cosθtilde1to2_for1

include("angle_functions.jl")
include("tbs_struct.jl")
include("rand_corr.jl")
include("rand_gen.jl")
include("lineshape.jl")
include("angular_functions.jl")
include("cross_channel_relations.jl")


export jls, twochain

export coupling_scheme12,
       coupling_scheme23,
       coupling_scheme31

export two_J,two_L,two_S,
       two_j,two_l,two_s

export posibleLS

include("coupling_scheme.jl")

end  # module ThreeBodyDecay
