module ThreeBodyDecay

using GSL
using QuadGK
# using WignerSymbols

export gσ1, gσ2, gσ3,
       σ3of1, σ1of2, σ2of3,
       cosθ12, cosθ23, cosθ31,
       cosθhat12,
       cosθhat23,
       cosθhat31,
       cos_plus_θhat31,
       cos_mins_θhat21,
       λ,
       Kibble, Kibble23, Kibble12, Kibble31

export ThreeBodySystem
export pthsq, pth

# density
export getbinned1dDensity,
       getbinned2dDensity
export flatDalitzPlotSample31
include("rand_corr.jl")
include("rand_gen.jl")

#lineshape
export pole, BW, BWdw
export Lineshape, BreitWigner, amp
export iRho, ChewMandestam
export iRhoQTB
include("lineshape.jl")

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

export cosζ31_for1,
       cosζ23_for1,
       cosζ12_for1

include("angle_functions.jl")
include("tbs_struct.jl")
include("angular_functions.jl")
include("cross_channel_relations.jl")

#
export jls, twochain

export coupling_scheme12,
       coupling_scheme23,
       coupling_scheme31

export two_J,two_L,two_S,
       two_j,two_l,two_s

export posibleLS

export HelicityRecoupling,
       HelicityRecoupling_doublearg
include("coupling_scheme.jl")

#
export amp_b2bzz,
    rate_b2bzz, rateCC_b2bzz, rateΛΛ_b2bzz,
    polSens_b2bzz
export get_couplings,
    set_couplings
include("half_to_half_zero_zero.jl")

end  # module ThreeBodyDecay
