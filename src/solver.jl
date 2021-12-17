function initialize!(qp::QP)
    idx = qp.idx

    @. qp.init.rhs[idx.x] = -qp.q
    @. qp.init.rhs[idx.s] = qp.h
    @. qp.init.rhs[qp.init.idx_y] = qp.b

    iterative_ref(qp.init.LS_F, qp.init.LS,
                           qp.init.sol,
                           qp.init.rhs,
                           qp.init.LS_IR)

    qp.x .= view(qp.init.sol,idx.x)
    qp.z .= view(qp.init.sol,idx.s)
    qp.y .= view(qp.init.sol,qp.init.idx_y)

    α = -Inf
    for i = 1:idx.ns
        if α < qp.z[i]
            α = qp.z[i]
        end
    end

    if α < 0
        for i = 1:idx.ns
            qp.s[i] = -qp.z[i]
        end
    else
        α += 1
        for i = 1:idx.ns
            qp.s[i] = -qp.z[i] + α
        end
    end

    α = -Inf
    for i = 1:idx.ns
        if α < -qp.z[i]
            α = -qp.z[i]
        end
    end

    if α >= 0
        α += 1
        for i = 1:idx.ns
            qp.z[i] = qp.z[i] + α
        end
    end


    return nothing
end
function solveqp!(qp::QP)

    @printf "iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n"
    @printf "---------------------------------------------------------\n"

    # @btime initialize!($qp)
    initialize!(qp)
    # @btime update_kkt!($qp)
    for i = 1:50

        # update all things KKT
        update_kkt!(qp)
        update_kkt_reg!(qp,1e-8)
        update_kkt_factor!(qp)

        # affine step
        rhs_kkt_a!(qp)
        iterative_ref(qp.KKT_F, qp.KKT, qp.p_a, qp.rhs_a, qp.KKT_IR)
        index_sol_a!(qp)
        # @show norm(qp.KKT\qp.rhs_a - qp.p_a)

        # centering and correcting step
        σ, μ = centering_params(qp)
        rhs_kkt_c!(qp, σ, μ)
        iterative_ref(qp.KKT_F, qp.KKT, qp.p_c, qp.rhs_c, qp.KKT_IR)
        index_sol_c!(qp)
        # @show norm(qp.KKT\qp.rhs_c - qp.p_c)
        # combine deltas
        combine_deltas!(qp)

        # last linesearch
        α = min(1,0.99*min(linesearch(qp.s,qp.Δ.s),linesearch(qp.z,qp.Δ.z)))

        update_vars!(qp,α)

        if logging(qp::QP,i,α)
            break
        end

    end
    return nothing
end
function quadprog(Q,q,A,b,G,h)
    qp = QP(Q,q,A,b,G,h)
    # @btime solveqp!($qp::QP)
    solveqp!(qp)
    return qp.x
end
