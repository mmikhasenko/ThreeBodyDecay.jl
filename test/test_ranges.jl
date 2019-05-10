
let m1 = 0.938, m2 = 0.49367, m3 = 0.13957, m0 = 2.46867
    tbs0 = ThreeBodySystem(m0,m1,m2,m3)
    @test tbs0.mthsq[1] == (m2+m3)^2
end
