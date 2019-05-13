# recoupling function
# two options for the type:
#   - either all arguments are Float64 or energy variables are complex
function change_basis_3from1(σ1, cosθ1, ϕ1, cosθ23, ϕ23,
        m1sq, m2sq, m3sq, s)
        # calculate σ3 in (23) frame
        m1 = sqrt(m1sq); sqrts = sqrt(s)
        sqrtλσ1 = sqrt((sqrts-m1)^2-σ1)*sqrt((sqrts+m1)^2-σ1)/(2*sqrt(σ1));
        σ3 = m1sq + m2sq +
        # 2*(  E2 * E1 -
            2*( (σ1 + m2sq - m3sq)/(2*sqrt(σ1)) * (s-m1sq-σ1)/(2*sqrt(σ1)) -
            # |p2|*cos(theta23) * (-|p1|)
            sqrt(λ(σ1, m2sq, m3sq)/(4*σ1))*cosθ23 * (-sqrtλσ1) );
        # caluclate p3 in (23) frame
        p23bu = sqrt(λ(σ1, m2sq, m3sq)/(4*σ1));
        p3_in23 = [(σ1 + m3sq - m2sq)/sqrt(4*σ1),
                   -p23bu*sqrt(1-cosθ23^2)*cos(ϕ23),
                   -p23bu*sqrt(1-cosθ23^2)*sin(ϕ23),
                   -p23bu*cosθ23];
        # boost to lab frame from (23) frame
        γ1 = (s + σ1 - m1sq)/sqrt(4*s*σ1);
        β1 = sqrt(1.0-1.0/γ1^2);
        p3_b = [γ1*(p3_in23[1]+β1*p3_in23[4]),
                p3_in23[2],
                p3_in23[3],
                γ1*(β1*p3_in23[1]+p3_in23[4])];

        ct1 = cosθ1; st1 = sqrt(1.0-cosθ1^2); cp1 = cos(ϕ1); sp1 = sin(ϕ1);
        # Rz(phi1) * Ry(theta1) * p3_boost
        p3_rot = [p3_b[1],
                  cp1*ct1* p3_b[2] + (-sp1)* p3_b[3] + cp1*st1* p3_b[4],
                  sp1*ct1* p3_b[2] +   cp1 * p3_b[3] + sp1*st1* p3_b[4],
                     -st1* p3_b[2] +     0 * p3_b[3] +     ct1* p3_b[4]];

        cosθ3 = -p3_rot[4]/sqrt(p3_b[2]^2+p3_b[3]^2+p3_b[4]^2);
        ϕ3 = (p3_rot[2] != zero(p3_rot[2])) ? atan(-p3_rot[3], -p3_rot[2]) : rand()*one(p3_rot[2]);

        cosθ12_n = m2sq + m3sq + 2* (σ3+m2sq-m1sq)/(2*sqrt(σ3)) * (s-m3sq-σ3)/(2*sqrt(σ3)) - σ1;
        m3 = sqrt(m3sq)
        sqrtλσ3 = sqrt((sqrts-m3)^2-σ3)*sqrt((sqrts+m3)^2-σ3)
        cosθ12_d = 2 * sqrt(λ(σ3, m1sq, m2sq)/(4*(σ3))) * sqrtλσ3 / (2*sqrt(σ3));
        cosθ12 = cosθ12_d ≈ 0.0+0.0im ? 2.0*rand()-1.00*one(σ1) : cosθ12_n/cosθ12_d;

        n1 = [0.,
             -sqrt(1.0-cosθ1^2)*cos(ϕ1),
             -sqrt(1.0-cosθ1^2)*sin(ϕ1),
             -cosθ1];
        ct3 = cosθ3; st3 = sqrt(1.0-cosθ3^2); cp3 = cos(ϕ3); sp3 = sin(ϕ3);
        # Rz(phi23) * Ry(theta23) * p3_boost
        n1_rot = [n1[1],
                  cp3*ct3* n1[2] + sp3*ct3* n1[3] + (-st3)* n1[4],
                   (-sp3)* n1[2] +    cp3 * n1[3] +   0.0 * n1[4],
                  cp3*st3* n1[2] + sp3*st3* n1[3] +   ct3 * n1[4]];
        # println(n1_rot[3], " ", n1_rot[2])
        ϕ12 = (n1_rot[2] != zero(n1_rot[2])) ? atan(n1_rot[3], n1_rot[2]) : rand()*one(n1_rot[2]);
        return σ3,cosθ3,ϕ3,cosθ12,ϕ12
end

change_basis_3from1(τ1, tbs::ThreeBodySystem) = change_basis_3from1(τ1..., tbs.msq[1],tbs.msq[2],tbs.msq[3], tbs.s)
change_basis_1from2(τ2, tbs::ThreeBodySystem) = change_basis_3from1(τ2..., tbs.msq[2],tbs.msq[3],tbs.msq[1], tbs.s)
change_basis_2from3(τ3, tbs::ThreeBodySystem) = change_basis_3from1(τ3..., tbs.msq[3],tbs.msq[1],tbs.msq[2], tbs.s)
