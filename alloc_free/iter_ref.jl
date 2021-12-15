using QDLDL, SparseArrays, SuiteSparse
using LinearAlgebra
using Infiltrator


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

# function ttttt()
#
#     n = 150
#     A = sprand(n,n,0.1)
#
#     A = A'*A + I
#
#     F = qdldl(A + 1e-4*I)
#     b = randn(n)
#     x = zeros(n)
#     r = zeros(n)
#     Δx = zeros(n)
#     @btime x = iterative_ref($F,$A,$x,$b,$r,$Δx)
#
#     @show norm(x - A\b)
#     return nothing
# end
#
# ttttt()
