using ThreeBodyDecay
using Test

ms = let
    mpi = 0.13957
    meta = 0.547862
    mp = 0.938
    # 
    ThreeBodyMasses(mpi, meta, mp; m0=19.0)
end

σs_b = border(ms)

@testset "Border is an array of MandestamTuple" begin
    @test σs_b isa Vector{T} where {T<:ThreeBodyDecay.MandestamTuple{Float64}}
    σs_b
end
