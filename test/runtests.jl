using SafeTestsets
using ThreeBodyDecay

# ThreeBodyDecay structure
@safetestset "Test of the tbs structure" begin
    include("test_threebodymasses.jl")
    include("test_threebodyspins.jl")
    include("test_threebodyparities.jl")
    include("test_invariants.jl")
end

@safetestset "Phases" begin include("test_phases.jl") end

include("wignerrotationdispatch.jl")

@safetestset "four vectors" begin include("test_generation.jl") end

# angular functions
@safetestset "Angular Functions block" begin
    include("circle_relations.jl")
    include("test_representaiton_property.jl")
    include("check_lineshape_functions.jl")
end

@safetestset "Permutations" begin include("test_wigner_permutations.jl") end

@safetestset "Sum rule" begin include("wigner_angle_sumrules.jl") end

# Lineshapes
@safetestset "Lineshapes" begin
    # include("breit-wigner.jl")
    include("scattering_length_X3872.jl")
end

# JPC functions
@safetestset "Couplings" begin include("test_couplings.jl") end

@safetestset "LS amplitudes" begin
    include("test_ls_amplitude.jl")
    include("test_sum_over_polarization.jl")
end

# Integrals
@safetestset "Integrals" begin
    include("test_three_body_phase_space_integral.jl")
end
# examples
@safetestset "Complete example" begin include("test_Xic_decay.jl") end

@safetestset "square dalitz" begin
    include("test_transformations_of_squaredalitz.jl")
end
