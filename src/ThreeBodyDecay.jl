module ThreeBodyDecay

using QuadGK
using StaticArrays
using Cuba
using PartialWaveFunctions
using Parameters
using RecipesBase
using Polynomials

export Kallen, sqrtKallenFact, Kibble,
        inphrange,
        lims, lims1, lims2, lims3
export phase
export ijk, ij_from_k
# 
export cosθij
export cosθ12, cosθ23, cosθ31
# 
export σkofi, σkofj
export σ1of2, σ2of3, σ3of1,
        σ1of3, σ2of1, σ3of2
include("kinematics.jl")

export x_str, minusone
include("utils.jl")

export ispositive
export TriavialWignerRotation
export wr
export cosζ
# 
export cosζ12_for0, cosζ23_for0, cosζ31_for0
# 
export cosζ21_for1, cosζ21_for2
export cosζ13_for1, cosζ13_for3
export cosζ32_for3, cosζ32_for2
# 
export cosζ12_for3, cosζ23_for1, cosζ31_for2
# 
export cosθhatk1, cosθhatk2, cosθhatk3
export cosζk1_for1, cosζk2_for2, cosζk3_for3
# 
include("wignerrotations.jl")


export ThreeBodyMasses, ThreeBodySpins, ThreeBodyParities
export MandestamTuple
export ThreeBodySystem
export DalitzPlotPoint
export Invariants
export randomPoint
export nt, x2, over2
export possible_helicities
export border31, border12, border23
export border13, border21, border32
export border
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
export jp, str2jp
export @jp_str
export ⊗
export possible_ls, possible_lsLS
export possible_coupling_schemes
include("coupling_scheme.jl")

export NoRecoupling, ParityRecoupling, RecouplingLS
export DecayChain, DecayChainLS, DecayChainsLS
export amplitude
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
