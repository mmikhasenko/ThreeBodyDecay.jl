using Test
using ThreeBodyDecay

using GSL
function ClGd_jsl(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    factor = sqrt(two_j+1)*(mod(two_j1-two_j2+two_m,4)==2 ? -1 : +1)
    three_j = sf_coupling_3j(two_j1,two_j2,two_j,two_m1,two_m2,-two_m)
    return factor*three_j;
end

for _ in 1:1_000_009
    two_j1 = rand(0:15); two_m1 = rand(-two_j1:2:two_j1)
    two_j2 = rand(0:15); two_m2 = rand(-two_j2:2:two_j2)
    two_j = rand(abs(two_j2-two_j1):2:(two_j1+two_j2))
    two_m = two_m1+two_m2
    v = ClGd_jsl(       two_j1,two_m1,two_j2,two_m2,two_j,two_m) -
        ClGd(two_j1,two_m1,two_j2,two_m2,two_j,two_m)
    (v > 0.1) && (@show (two_j1,two_m1,two_j2,two_m2,two_j,two_m))
end
