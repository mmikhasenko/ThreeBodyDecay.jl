using SafeTestsets
using ThreeBodyDecay

# ThreeBodyDecay structure
@safetestset "Test of the tbs structure" begin
    include("threebodymasses.jl")
    include("threebodyspins.jl")
    include("threebodyparities.jl")
    include("invariants.jl")
end

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

# Lineshapes
@safetestset "Lineshapes" begin
    # include("breit-wigner.jl")
    include("scattering_length_X3872.jl")
end

# JPC functions
@safetestset "Couplings" begin include("test_couplings.jl") end

@safetestset "LS amplitudes" begin
    include("ls_amplitude.jl")
    include("sum_over_polarization.jl")
end

# Integrals
@safetestset "Integrals" begin
    include("three_body_phase_space_integral.jl")
end
# examples
@safetestset "Complete example" begin include("test_Xic_decay.jl") end

@safetestset "square dalitz" begin
    include("transformations_of_squaredalitz.jl")
end
