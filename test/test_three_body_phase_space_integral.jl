using ThreeBodyDecay
using Test
using Parameters

function QTB_mismatch_factor(dc)
  k = dc.k
  i, j = ij_from_k(k)
  @unpack Hij, HRk = dc
  @unpack two_s, tbs = dc
  @unpack two_js = tbs
  #
  avHHsq =
    sum((tbs.two_js[4] != (two_τ - two_λs[k]) ? 0.0 : 1.0) *
        abs2(amplitude(Hij, two_τ, two_λs[k])) *
        abs2(amplitude(Hij, two_λs[i], two_λs[j]))
        for two_λs in itr(tbs.two_js),
        two_τ in -dc.two_s:2:dc.two_s)
  return avHHsq
end

@testset "QaudGK vs Cuba within 0.01%" begin
  ms = ThreeBodyMasses(m0=6.0, m1=1.0, m2=1.5, m3=2.0)
  #
  v_qtv = RhoQTB(ms.m0^2, σ -> 1.0, ms^2)
  v_fll = real(three_body_phase_space_integral(σs -> 1.0, ms))
  #
  @test (v_qtv - v_fll) / v_fll < 1e-5
end

@testset "QaudGK vs Cuba for P-wave within 0.01%" begin
  tbs = ThreeBodySystem(1.0, 1.5, 2.0, m0=6.0)
  #
  ch = DecayChainLS(1, σ -> 1.0; two_s=2, tbs=tbs, parity='-', Ps=['+', '+', '+', '+'])
  sp = summed_over_polarization((σs, two_λ) -> abs2(amplitude(ch, σs, two_λ)), tbs.two_js)
  #
  v_qtv = RhoQTB(tbs.ms.m0^2, σ -> 1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  #
  @test abs(v_qtv - v_fll) / v_fll < 1e-5
end


tbs = ThreeBodySystem(1.0, 1.5, 2.0; m0=6.0, two_js=ThreeBodySpins(0, 1, 0; two_h0=1))


@testset "QaudGK vs Cuba for half-integer spin within 0.01%" begin
  #
  ch = DecayChainLS(1, σ -> 1.0; two_s=1, tbs=tbs, parity='-', Ps=['-', '+', '-', '+'])
  sp = summed_over_polarization((σs, two_λ) -> abs2(amplitude(ch, σs, two_λ)), tbs.two_js)
  #
  v_qtv = RhoQTB(tbs.ms.m0^2, σ -> 1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  # @show v_fll
  #
  @test abs(v_qtv - v_fll) / v_fll < 1e-5
end

@testset "QaudGK vs Cuba for half-integer spin and P-wave within 0.01%" begin
  #
  ch = DecayChainLS(1, σ -> 1.0; two_s=1, tbs=tbs, parity='+', Ps=['-', '+', '-', '+'])
  sp = summed_over_polarization((σs, two_λ) -> abs2(amplitude(ch, σs, two_λ)), tbs.two_js)
  # sp = σs->abs2(amplitude(ch,σs,[0,0,0,0]))
  #
  v_qtv = RhoQTB(tbs.ms.m0^2, σ -> 1.0, tbs.ms^2)
  v_qtv *= QTB_mismatch_factor(ch) # helicity couplings
  # @show v_qtv
  v_fll = real(three_body_phase_space_integral(sp, tbs.ms))
  # @show v_fll
  #
  @test abs(v_qtv - v_fll) / v_fll < 1e-5
end
