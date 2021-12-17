function rhs_kkt_a!(qp::QP)
    idx = qp.idx
    c = qp.cache

    # qp.rhs_a[idx.x] = -(qp.A'*qp.y + qp.G'*qp.z + qp.Q*qp.x + qp.q)
    mul!(c.x.c1,qp.A',qp.y)
    mul!(c.x.c2,qp.G',qp.z)
    mul!(c.x.c3,qp.Q,qp.x)
    @. qp.rhs_a[idx.x] = -c.x.c1 - c.x.c2 - c.x.c3 - qp.q

    # qp.rhs_a[idx.s] = -(qp.z)
    @. qp.rhs_a[idx.s] = -qp.z

    # qp.rhs_a[idx.z] = -(qp.G*qp.x + qp.s - qp.h)
    mul!(c.z.c1,qp.G,qp.x)
    @. qp.rhs_a[idx.z] = -c.z.c1 - qp.s + qp.h

    # qp.rhs_a[idx.y] = -(qp.A*qp.x - qp.b)
    mul!(c.y.c1,qp.A,qp.x)
    @. qp.rhs_a[idx.y] = -c.y.c1 + qp.b
    return nothing
end
function index_sol_a!(qp::QP)
    # qp.Δ.x_a .= qp.p_a[qp.idx.x]
    qp.Δ.x_a .= view(qp.p_a,qp.idx.x)

    # qp.Δ.s_a .= qp.p_a[qp.idx.s]
    qp.Δ.s_a .= view(qp.p_a,qp.idx.s)

    # qp.Δ.z_a .= qp.p_a[qp.idx.z]
    qp.Δ.z_a .= view(qp.p_a,qp.idx.z)

    # qp.Δ.y_a .= qp.p_a[qp.idx.y]
    qp.Δ.y_a .= view(qp.p_a,qp.idx.y)
    return nothing
end
function index_sol_c!(qp::QP)
    qp.Δ.x_c .= view(qp.p_c,qp.idx.x)
    qp.Δ.s_c .= view(qp.p_c,qp.idx.s)
    qp.Δ.z_c .= view(qp.p_c,qp.idx.z)
    qp.Δ.y_c .= view(qp.p_c,qp.idx.y)
    return nothing
end

function rhs_kkt_c!(qp::QP, σ, μ)
    idx = qp.idx
    qp.rhs_c .= 0
    # qp.rhs_c[idx.s] = (σ*μ .- (qp.Δ.s_a .* qp.Δ.z_a)) ./ qp.s
    @. qp.rhs_c[idx.s] = (σ*μ - (qp.Δ.s_a * qp.Δ.z_a)) / qp.s
    return nothing
end

function initialize_kkt!(KKT,idx::IDX,Q,G,A)
    KKT[idx.x, idx.x] = Q
    KKT[idx.x, idx.z] = G'
    KKT[idx.x, idx.y] = A'
    KKT[idx.s, idx.s] = Diagonal(ones(idx.nz))
    KKT[idx.s, idx.z] = Diagonal(ones(idx.ns))
    KKT[idx.z, idx.x] = G
    KKT[idx.z, idx.s] = I(idx.nz)
    KKT[idx.y, idx.x] = A
    return nothing
end
function regularize!(A,np,nd,ρ)
    for i = 1:np
        A[i,i] += ρ
    end
    for i = (np + 1):(np + nd)
        A[i,i] -= ρ
    end
    return nothing
end

@inline function update_kkt_factor!(qp::QP)
    QDLDL.update_A!(qp.KKT_F,qp.KKT_reg)
    return nothing
end
function update_kkt_reg!(qp::QP,ρ)
    idx = qp.idx
    # qp.KKT[idx.s, idx.s] = Diagonal(qp.z ./ qp.s)
    # qp.KKT[idx.s, idx.z] = I(idx.ns)
    for i = 1:idx.ns
        qp.KKT_reg[idx.s[i],idx.s[i]] = qp.KKT[idx.s[i],idx.s[i]] + ρ
    end
    return nothing
end

function update_kkt!(qp::QP)
    idx = qp.idx
    # qp.KKT[idx.s, idx.s] = Diagonal(qp.z ./ qp.s)
    # qp.KKT[idx.s, idx.z] = I(idx.ns)
    for i = 1:idx.ns
        qp.KKT[idx.s[i],idx.s[i]] = qp.z[i]/qp.s[i]
        qp.KKT[idx.s[i],idx.z[i]] = 1.0
    end
    return nothing
end
