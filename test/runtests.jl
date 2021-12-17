using QuadProg
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Test

include("create_MPC.jl")


@testset "mpc test" begin
    N = 15
    Q,q,A,b,G,h = create_MPC(N)
    x1 = quadprog(Q,q,A,b,G,h)
end
