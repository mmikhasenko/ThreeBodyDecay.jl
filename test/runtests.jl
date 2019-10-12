using ThreeBodyDecay
using Test

# ThreeBodyDecay structure
include("test_ranges.jl")

# angular functions
include("circle_relations.jl")
include("representaiton_property.jl")
include("check_lineshape_functions.jl")

# JPC functions
include("test_couplings.jl")

# examples
include("test_Xic_decay.jl")
