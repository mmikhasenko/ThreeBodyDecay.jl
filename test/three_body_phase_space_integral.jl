using Test
using ThreeBodyDecay

@testset "QaudGK vs Cuba within 0.01%" begin
  tbs = ThreeBodySystem(6,1,1.5,2)
  #
  v_qtv = RhoQTB(tbs.msq[4],(s,σ)->1.0, tbs.msq)
  v_fll = real(three_body_phase_space_integral(σs->1.0, tbs))
  #
  @test (v_qtv-v_fll)/v_fll < 1e-5
end

@testset "QaudGK vs Cuba for P-wave within 0.01%" begin
  tbs = ThreeBodySystem(6,1,1.5,2, two_jps=([0,0,0,0],['+','+','+','+']))
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=2, tbs=tbs, parity='-')
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs)
  #
  v_qtv = RhoQTB(tbs.msq[4],(s,σ)->1.0, tbs.msq)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs))
  @show v_fll
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end


@testset "QaudGK vs Cuba for half-integer spin within 0.01%" begin
  tbs = ThreeBodySystem(6,1,1.5,2, two_jps=([0, 1, 0, 1], ['-','+','-','+']))
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=1, tbs=tbs, parity='-')
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs)
  #
  v_qtv = RhoQTB(tbs.msq[4],(s,σ)->1.0, tbs.msq)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs))
  @show v_fll
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end

@testset "QaudGK vs Cuba for half-integer spin and P-wave within 0.01%" begin
  tbs = ThreeBodySystem(6,1,1.5,2, two_jps=([0, 1, 0, 1], ['-','+','-','+']))
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=1, tbs=tbs, parity='+')
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs)
  # sp = σs->abs2(amplitude(σs,[0,0,0,0],ch))
  #
  v_qtv = RhoQTB(tbs.msq[4],(s,σ)->1.0, tbs.msq)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs))
  @show v_fll
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end
