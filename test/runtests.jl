using ThreeBodyDecay
using Test

@testset "ThreeBodyDecay.jl" begin

    include("test_ranges.jl")
    include("circle_relations.jl")
    include("representaiton_property.jl")
    include("check_lineshape_functions.jl")

end
