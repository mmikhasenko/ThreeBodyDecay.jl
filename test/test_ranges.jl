
const m1 = 0.938
const m2 = 0.49367
const m3 = 0.13957
const m0 = 2.46867

tbs0 = ThreeBodySystem(m0,m1,m2,m3)

@test tbs0.mthsq[1] == (m2+m3)^2
