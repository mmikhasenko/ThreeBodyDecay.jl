using SafeTestsets

# ThreeBodyDecay structure
@safetestset "Test of the tbs structure" begin include("test_ranges.jl") end

@safetestset "four vectors" begin include("test_generation.jl") end

# angular functions
@safetestset "Angular Functions block" begin
    include("circle_relations.jl")
    include("representaiton_property.jl")
    include("check_lineshape_functions.jl")
end

@safetestset "Permutations" begin include("wigner_permutations.jl") end
@safetestset "Phases" begin include("phases.jl") end

@safetestset "Sum rule" begin include("wigner_angle_sumrules.jl") end

# JPC functions
@safetestset "Couplings" begin include("test_couplings.jl") end

# examples
@safetestset "Complete example" begin include("test_Xic_decay.jl") end
