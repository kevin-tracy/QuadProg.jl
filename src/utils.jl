function combine_deltas!(qp::QP)
    @. qp.Δ.x = qp.Δ.x_a + qp.Δ.x_c
    @. qp.Δ.s = qp.Δ.s_a + qp.Δ.s_c
    @. qp.Δ.z = qp.Δ.z_a + qp.Δ.z_c
    @. qp.Δ.y = qp.Δ.y_a + qp.Δ.y_c
    return nothing
end
function update_vars!(qp::QP,α::Float64)
    @. qp.x += α*qp.Δ.x
    @. qp.s += α*qp.Δ.s
    @. qp.z += α*qp.Δ.z
    @. qp.y += α*qp.Δ.y
    return nothing
end

function logging(qp::QP,iter,α)

    c = qp.cache

    # J = 0.5*qp.x'*qp.Q*qp.x + dot(qp.q,qp.x)
    c.x.c1 .= 0
    mul!(c.x.c1,qp.Q,qp.x)
    J = 0.5*dot(qp.x,c.x.c1) + dot(qp.q,qp.x)

    # duality gap
    gap = dot(qp.s,qp.z)/qp.idx.ns

    # |Ax - b|
    c.y.c1 .= 0
    mul!(c.y.c1,qp.A,qp.x)
    @. c.y.c1 -= qp.b
    eq_res = norm(c.y.c1)

    # |Gx + s - h|
    c.z.c1 .= 0
    mul!(c.z.c1,qp.G,qp.x)
    @. c.z.c1 += qp.s - qp.h
    ineq_res = norm(c.z.c1)

    @printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter, J, gap, eq_res,
          ineq_res, α)

    return (gap<1e-8)
end
