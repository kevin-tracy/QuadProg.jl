using LinearAlgebra, SparseArrays
using SuiteSparse
# using Infiltrator
# using QDLDL, ECOS, Convex
using QDLDL
using Convex, Mosek, MosekTools, Printf
# using Convex, ECOS
using Printf
using Infiltrator
using Test

include(joinpath(@__DIR__,"create_MPC.jl"))

struct x_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function x_cache(nx)
        new(zeros(nx),zeros(nx),zeros(nx))
    end
end
struct z_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function z_cache(nz)
        new(zeros(nz),zeros(nz),zeros(nz))
    end
end
struct y_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function y_cache(ny)
        new(zeros(ny),zeros(ny),zeros(ny))
    end
end
struct CACHE
    x::x_cache
    z::z_cache
    y::y_cache
    function CACHE(nx,nz,ny)
        new(x_cache(nx),z_cache(nz),y_cache(ny))
    end
end
struct IDX
    # contains the variable indexing stuff for [x;s;z;y]
    nx::Int64
    ns::Int64
    nz::Int64
    ny::Int64
    N::Int64

    x::UnitRange{Int64}
    s::UnitRange{Int64}
    z::UnitRange{Int64}
    y::UnitRange{Int64}
end
struct IterRef
    r::Vector{Float64}
    Δx::Vector{Float64}
    function IterRef(N)
        new(zeros(N),zeros(N))
    end
end
struct INIT
    LS::SparseMatrixCSC{Float64, Int64}
    LS_reg::SparseMatrixCSC{Float64, Int64}
    LS_F::QDLDL.QDLDLFactorisation{Float64, Int64}
    rhs::Vector{Float64}
    sol::Vector{Float64}
    idx_y::UnitRange{Int64}
    LS_IR::IterRef
    function INIT(nx,nz,ny,idx,Q,G,A,idx_y)
        Ni = nx + nz + ny
        LS = spzeros(Ni, Ni)
        LS[idx.x,idx.x] = Q
        LS[idx.x,idx.s] = G'
        LS[idx.x,idx_y] = A'
        LS[idx.s,idx.x] = G
        LS[idx.s,idx.s] = -I(idx.ns)
        LS[idx_y,idx.x] = A
        LS_reg = copy(LS)
        regularize!(LS_reg, nx, nz + ny,1e-8)
        new(LS, LS_reg, qdldl(LS_reg), zeros(Ni), zeros(Ni),idx_y,IterRef(Ni))
    end
end
struct DELTA
    # stores all the deltas (search directions)

    # affine deltas
    x_a::Array{Float64,1}
    s_a::Array{Float64,1}
    z_a::Array{Float64,1}
    y_a::Array{Float64,1}

    # centering + correcting deltas
    x_c::Array{Float64,1}
    s_c::Array{Float64,1}
    z_c::Array{Float64,1}
    y_c::Array{Float64,1}

    # total deltas
    x::Array{Float64,1}
    s::Array{Float64,1}
    z::Array{Float64,1}
    y::Array{Float64,1}

    # constructor
    function DELTA(nx,ns,nz,ny)
        new(zeros(nx),zeros(ns),zeros(nz),zeros(ny),
            zeros(nx),zeros(ns),zeros(nz),zeros(ny),
            zeros(nx),zeros(ns),zeros(nz),zeros(ny))
    end
end

struct QP

    # problem data
    Q::SparseMatrixCSC{Float64,Int64}
    q::Array{Float64,1}
    A::SparseMatrixCSC{Float64,Int64}
    b::Array{Float64,1}
    G::SparseMatrixCSC{Float64,Int64}
    h::Array{Float64,1}

    # variables
    x::Array{Float64,1} # primal
    s::Array{Float64,1} # primal
    z::Array{Float64,1} # dual
    y::Array{Float64,1} # dual

    # KKT stuff
    KKT::SparseMatrixCSC{Float64,Int64}
    KKT_reg::SparseMatrixCSC{Float64,Int64}
    KKT_F::QDLDL.QDLDLFactorisation{Float64, Int64}
    KKT_IR::IterRef
    rhs_a::Array{Float64,1}
    rhs_c::Array{Float64,1}
    p_a::Array{Float64,1}
    p_c::Array{Float64,1}

    # indexing
    idx::IDX

    # deltas
    Δ::DELTA

    # cache
    cache::CACHE

    # initialization stuff
    init::INIT

    function QP(Q,q,A,b,G,h)

        # length of variables
        nx = length(q)
        ns = length(h)
        nz = length(h)
        ny = length(b)
        N = nx + ns + nz + ny

        # indexing when stacked [x;s;z;y]
        idx_x = 1:nx
        idx_s = (nx + 1) : (nx + ns)
        idx_z = (nx + ns + 1) : (nx + ns + nz)
        idx_y = (nx + ns + nz + 1) : (nx + ns + nz + ny)
        idx = IDX(nx,ns,nz,ny,N, idx_x, idx_s, idx_z, idx_y)

        # kkt stuff
        KKT = spzeros(N,N)
        initialize_kkt!(KKT,idx,Q,G,A)
        KKT_reg = copy(KKT)
        regularize!(KKT_reg, nx + ns,nz + ny, 1e-8) # TODO: update this ρ
        KKT_F = qdldl(KKT_reg)
        KKT_IR = IterRef(N)
        rhs_a = zeros(N)
        rhs_c = zeros(N)
        p_a = zeros(N)
        p_c = zeros(N)

        # initialize variables to zero
        x = zeros(nx)
        s = zeros(ns)
        z = zeros(nz)
        y = zeros(ny)

        # deltas
        Δ = DELTA(nx,ns,nz,ny)

        #initialization
        init = INIT(nx,nz,ny,idx,Q,G,A,idx_y .- ns)
        # nx,nz,ny,idx,Q,G,A,idx_y

        new(Q,q,A,b,G,h,x,s,z,y, KKT, KKT_reg, KKT_F, KKT_IR, rhs_a, rhs_c, p_a, p_c, idx, Δ,CACHE(nx,nz,ny),init)
    end


end

# ---------------real functions---------------
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
function linesearch(x,dx)
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

function rhs_kkt_c!(qp::QP, σ, μ)
    idx = qp.idx
    qp.rhs_c .= 0
    # qp.rhs_c[idx.s] = (σ*μ .- (qp.Δ.s_a .* qp.Δ.z_a)) ./ qp.s
    @. qp.rhs_c[idx.s] = (σ*μ - (qp.Δ.s_a * qp.Δ.z_a)) / qp.s
    return nothing
end
function combine_deltas!(qp::QP)
    @. qp.Δ.x = qp.Δ.x_a + qp.Δ.x_c
    @. qp.Δ.s = qp.Δ.s_a + qp.Δ.s_c
    @. qp.Δ.z = qp.Δ.z_a + qp.Δ.z_c
    @. qp.Δ.y = qp.Δ.y_a + qp.Δ.y_c
    return nothing
end
function update_vars!(qp::QP,α)
    @. qp.x += α*qp.Δ.x
    @. qp.s += α*qp.Δ.s
    @. qp.z += α*qp.Δ.z
    @. qp.y += α*qp.Δ.y
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

    # # initialize_kkt!(qp)
    #
    # for i = 1:50
    #
    #     # update linear system for solves
    #     update_kkt!(qp)
    #     # kkt_factor = qdldl(qp.KKT +
    #                       # 1e-10*Diagonal([ones(qp.idx.nx + qp.idx.ns);-ones(qp.idx.nz + qp.idx.ny)]))
    #
    #     # mat"figure
    #     # hold on
    #     # spy($qp.KKT)
    #     # hold off"
    #     # @show size(qp.KKT)
    #     # @show rank(qp.KKT)
    #     # @show cond(Array(qp.KKT),2)
    #     # error()
    #     # @infiltrate
    #     kkt_factor = lu(qp.KKT)
    #
    #     # affine step
    #     rhs_kkt_a!(qp)
    #
    #
    #     qp.p_a .= kkt_factor\qp.rhs_a
    #     # qp.p_a .= iterative_ref(kkt_factor,qp.KKT,qp.rhs_a)
    #     index_sol_a!(qp)
    #
    #
    #     # centering and correcting step
    #     σ, μ = centering_params(qp)
    #
    #     rhs_kkt_c!(qp, σ, μ)
    #     qp.p_c .= kkt_factor\qp.rhs_c
    #     # qp.p_c .= iterative_ref(kkt_factor,qp.KKT,qp.rhs_c)
    #     index_sol_c!(qp)
    #
    #     # combine deltas
    #     combine_deltas!(qp)
    #
    #     # last linesearch
    #     α = min(1,0.99*min(linesearch(qp.s,qp.Δ.s),linesearch(qp.z,qp.Δ.z)))
    #
    #     update_vars!(qp,α)
    #
    #     if logging(qp::QP,i,α)
    #         break
    #     end
    # end
    # if dot(qp.s,qp.z) > 1e-8
    #     error("PDIP FAILED")
    # end
    return nothing
end

function logging(qp::QP,iter,α)

    J = 0.5*qp.x'*qp.Q*qp.x + dot(qp.q,qp.x)
    gap = dot(qp.s,qp.z)
    eq_res = norm(qp.A*qp.x - qp.b)
    ineq_res = norm(qp.G*qp.x + qp.s - qp.h)


    @printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter, J, gap, eq_res,
          ineq_res, α)

    return (gap<1e-8)
end

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

# function initialize_kkt!(qp::QP)
#     idx = qp.idx
#
#     qp.KKT[idx.x, idx.x] = qp.Q
#     qp.KKT[idx.x, idx.z] = qp.G'
#     qp.KKT[idx.x, idx.y] = qp.A'
#     qp.KKT[idx.s, idx.s] = Diagonal(qp.z)
#     qp.KKT[idx.s, idx.z] = Diagonal(qp.s)
#     qp.KKT[idx.z, idx.x] = qp.G
#     qp.KKT[idx.z, idx.s] = I(idx.nz)
#     qp.KKT[idx.y, idx.x] = qp.A
#
#     return nothing
# end
function initialize_kkt!(KKT,idx::IDX,Q,G,A)
    # idx = qp.idx

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

function quadprog(Q,q,A,b,G,h)
    qp = QP(Q,q,A,b,G,h)
    # @btime solveqp!($qp::QP)
    solveqp!(qp)
    return qp.x
end

function tt()
    N = 20
    Q,q,A,b,G,h = create_MPC(N)
    x1 = quadprog(Q,q,A,b,G,h)
    x = Variable(length(q))
    problem = minimize(0.5*quadform(x,Matrix(Q)) + dot(q,x),[A*x == b, G*x <= h])
    Convex.solve!(problem,Mosek.Optimizer)
    @test norm(x.value - x1) <1e-3
    #
    # # @show length(x1)
    # X = reshape(x1[1:end-6],9,N-1)
    # U = X[7:9,:]
    # mat"
    # figure
    # hold on
    # plot($U')
    # hold off
    # "
end


# create_MPC()
tt()