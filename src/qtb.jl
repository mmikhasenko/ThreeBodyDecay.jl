# 
function QTB_mismatch_factor(dc)
    k = dc.k; i,j = ij_from_k(k);
    tbs = dc.tbs
    two_js = tbs.two_js
    two_s = dc.two_s
    two_J = dc.tbs.two_js[0]
    #
    avHHsq =
        sum((tbs.two_js[4]!=(two_τ-two_λs[k]) ? 0.0 : 1.0) *
            jls_coupling(two_js[i], two_λs[i], two_js[j], two_λs[j], two_s, dc.two_ls[1], dc.two_ls[2])^2 *
            jls_coupling(two_s, two_τ, two_js[k], two_λs[k], two_js[4], dc.two_LS[1], dc.two_LS[2])^2
            for two_λs in itr(tbs.two_js),
                two_τ in -dc.two_s:2:dc.two_s)
    return (two_J+1) * avHHsq
end
