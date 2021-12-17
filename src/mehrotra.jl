function linesearch(x::Vector{Float64},dx::Vector{Float64})
    # α = min(1.0, minimum([dx[i]<0 ? -x[i]/dx[i] : Inf for i = 1:length(x)]))
    α = 1.0
    for i = 1:length(x)
        if dx[i]<0
            α = min(α,-x[i]/dx[i])
        end
    end
    return α
end
function centering_params(qp::QP)

    μ = dot(qp.s,qp.z)/qp.idx.ns

    α = min(linesearch(qp.s,qp.Δ.s_a), linesearch(qp.z,qp.Δ.z_a))

    # σ = (dot(qp.s + α*qp.Δ.s_a, qp.z + α*qp.Δ.z_a)/dot(qp.s,qp.z))^3
    @. qp.cache.z.c1 = qp.s + α*qp.Δ.s_a
    @. qp.cache.z.c2 = qp.z + α*qp.Δ.z_a
    σ = (dot(qp.cache.z.c1,qp.cache.z.c2)/dot(qp.s,qp.z))^3

    return σ, μ
end
