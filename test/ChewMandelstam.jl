using Test
using ThreeBodyDecay

let m1sq = 1.0, m2sq = 2.0, iϵ = 1e-8im
  sv = (√m1sq+√m2sq)^2 .+ rand(10) # rand -> [0, 1]
  v1 = iRho.(sv .+ iϵ, 1.0, 2.0)
  v2 = 1im .* imag.(ChewMandestam.(sv .+ iϵ,1.0,2.0))
  @test max(abs.(v1 - v2)...) < 10*imag(iϵ)
end
