struct IterRef
    r::Vector{Float64}
    Δx::Vector{Float64}
    function IterRef(N)
        new(zeros(N),zeros(N))
    end
end
function iterative_ref(F::QDLDL.QDLDLFactorisation{Float64, Int64},
                       A::SparseMatrixCSC{Float64, Int64},
                       x::Vector{Float64},
                       b::Vector{Float64},
                       IR::IterRef)
    x .= 0
    @. IR.r = b
    for i = 1:10
        @. IR.Δx = IR.r
        QDLDL.solve!(F,IR.Δx)
        @. x += IR.Δx

        # this is me using Δx as a lazy cache
        mul!(IR.Δx,A,x)
        @. IR.r = b - IR.Δx
        # if dot(IR.r,IR.r)<1e-25
        #     break
        # end
    end
    return nothing
end
