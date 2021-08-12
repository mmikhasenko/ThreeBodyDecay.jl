using Test
using ThreeBodyDecay

@testset "QaudGK vs Cuba within 0.01%" begin
  ms = ThreeBodyMasses(m0=6,m1=1,m2=1.5,m3=2)
  #
  v_qtv = RhoQTB(ms.m0^2,(s,σ)->1.0, ms^2)
  v_fll = real(three_body_phase_space_integral(σs->1.0, ms))
  #
  @test (v_qtv-v_fll)/v_fll < 1e-5
end

@testset "QaudGK vs Cuba for P-wave within 0.01%" begin
  tbs = ThreeBodySystem(1,1.5,2, m0=6)
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=2, tbs=tbs, parity='-', Ps=['+','+','+','+'])
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs.two_js)
  #
  v_qtv = RhoQTB(tbs.ms.m0^2,(s,σ)->1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end


tbs = ThreeBodySystem(1,1.5,2; m0=6, two_js=ThreeBodySpins(0, 1, 0; two_h0=1))

@testset "QaudGK vs Cuba for half-integer spin within 0.01%" begin
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=1, tbs=tbs, parity='-', Ps=['-','+','-','+'])
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs.two_js)
  #
  v_qtv = RhoQTB(tbs.ms.m0^2,(s,σ)->1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  # @show v_fll
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end

@testset "QaudGK vs Cuba for half-integer spin and P-wave within 0.01%" begin
  #
  ch = decay_chain(1,(s,σ)->1.0; two_s=1, tbs=tbs, parity='+', Ps=['-','+','-','+'])
  sp = summed_over_polarization((σs,two_λ)->abs2(amplitude(σs,two_λ,ch)), tbs.two_js)
  # sp = σs->abs2(amplitude(σs,[0,0,0,0],ch))
  #
  v_qtv = RhoQTB(tbs.ms.m0^2,(s,σ)->1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  # @show v_fll
  #
  @test abs(v_qtv-v_fll)/v_fll < 1e-5
end
