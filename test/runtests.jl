using ThreeBodyDecay
using Test

@testset "ThreeBodyDecay.jl" begin

    include("test_ranges.jl")
    include("circle_relations.jl")
    include("representaiton_property.jl")

end
