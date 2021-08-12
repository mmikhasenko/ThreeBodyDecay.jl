module ThreeBodyDecay

using QuadGK
using StaticArrays
using Cuba
using PartialWaveFunctions
using Parameters
using RecipesBase

export gσ1, gσ2, gσ3,
       σkofi, σ3of1, σ1of2, σ2of3,
       cosθ12, cosθ23, cosθ31,
       cosθhat12,
       cosθhat23,
       cosθhat31,
       cos_plus_θhat31,
       cos_mins_θhat21,
       λ,
       Kibble,
       inphrange,
       lims,lims1,lims2,lims3
export cosζ13_for1,
       cosζ21_for1
export cosζ23_for1, cosζ31_for2, cosζ12_for3
export cosζk1_for1,
       cosζk2_for2,
       cosζk3_for3,
       cosθhatk1,
       cosθij
export phase
export ij_from_k
include("angle_functions.jl")

export ThreeBodyMasses, ThreeBodySpins, ThreeBodyParities
export ThreeBodySystem
export DalitzPlotPoint
export Invariants
export randomPoint
export nt, x2, over2
export possible_helicities
export border, border31, border12, border23
export flatDalitzPlotSample
include("tbs_struct.jl")

# density
export getbinned1dDensity,
       getbinned2dDensity
export gridded_density_function
include("rand_corr.jl")

export rotz!, roty!, roty_cos!, roty_cos_inv!, boostz!
export invmasssq
export four_vectors_in_binary_decay
include("rand_gen.jl")

#lineshape
export pole, BW, BWdw
export Lineshape, BreitWigner, amp
export ScattLenApproximation, ScattLen
export Rho, iRho, ChewMandestam
export RhoQTB, iRhoQTB
#
export three_body_phase_space_integral
include("lineshape.jl")

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

export HelicityRecoupling,
       HelicityRecoupling_doublearg
export jp, @jp_str
export ⊗
export possible_ls, possible_lsLS
export possible_coupling_schemes
include("coupling_scheme.jl")

export Zksτ
include("general_case.jl")

export QTB_mismatch_factor
include("qtb.jl")

export decay_chain, decay_chains, amplitude
export itr, summed_over_polarization
include("decay_channel.jl")

include("dalitzplotsrecipe.jl")

export squaredalitz,
       squaredalitz1,
       squaredalitz2,
       squaredalitz3
export invsquaredalitz,
       invsquaredalitz1,
       invsquaredalitz2,
       invsquaredalitz3
export jacobean_squaredalitz,
       jacobean_squaredalitz1,
       jacobean_squaredalitz2,
       jacobean_squaredalitz3
include("squaredalitz.jl")

end  # module ThreeBodyDecay
