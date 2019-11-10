module ThreeBodyDecay

using GSL
using QuadGK
using StaticArrays
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
export cosζ13_for1,
       cosζ23_for1,
       cosζ21_for1
export cosζk1_for1,
       cosζk2_for2,
       cosζk3_for3,
       cosθhatk1,
       cosθij
export phase
export ij_from_k
include("angle_functions.jl")


export ThreeBodySystem, DalitzPlotPoint,
       DalitzPlotPoint12, DalitzPlotPoint23, DalitzPlotPoint31,
       randomPoint
export s,m0,m1,m2,m3
export two_J, two_Λ
export possible_helicities
include("tbs_struct.jl")

# density
export getbinned1dDensity,
       getbinned2dDensity
export flatDalitzPlotSample31
include("rand_corr.jl")

export rotz!, roty!, roty_cos!, boostz!
export invmasssq
export four_vectors_in_binary_decay
include("rand_gen.jl")

#lineshape
export pole, BW, BWdw
export Lineshape, BreitWigner, amp
export Rho, iRho, ChewMandestam
export RhoQTB, iRhoQTB
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
include("angular_functions.jl")

export change_basis_3from1,
        change_basis_1from2,
        change_basis_2from3

include("cross_channel_relations.jl")

#
export jls, twochain

export coupling_scheme, coupling_schemek,
       coupling_scheme12, coupling_scheme23, coupling_scheme31

export clebsch_for_chaink, jls_coupling

export two_J,two_L,two_S,
       two_j,two_l,two_s

export possibleLS

export HelicityRecoupling,
       HelicityRecoupling_doublearg
include("coupling_scheme.jl")

#
export Zksτ
export decay_chain, amplitude
export summed_over_polarization
include("general_case.jl")

#
export amp_b2bzz,
    rate_b2bzz, rateCC_b2bzz, rateΛΛ_b2bzz,
    polSens_b2bzz
export get_couplings,
    set_couplings
include("half_to_half_zero_zero.jl")

end  # module ThreeBodyDecay
