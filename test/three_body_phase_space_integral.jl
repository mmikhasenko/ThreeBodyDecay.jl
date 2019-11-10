using Test
using ThreeBodyDecay

@testset "QaudGK integrator vs Cuba within 0.01%" begin
  tbs = ThreeBodySystem(6,1,1.5,2)
  #
  v_qtv = RhoQTB(tbs.msq[4],(s,σ)->1.0, tbs.msq)
  v_fll = real(three_body_phase_space_integral(σs->1.0, tbs))
  #
  @test (v_qtv-v_fll)/v_fll < 1e-5
end
