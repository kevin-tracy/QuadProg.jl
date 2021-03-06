using LinearAlgebra, SparseArrays
using SuiteSparse
# using Infiltrator
# using QDLDL, ECOS, Convex
using QDLDL
# using Convex, Mosek, MosekTools, Printf
using OSQP
# using Convex, ECOS
using Printf
using Infiltrator
using Test
using BenchmarkTools

include(joinpath(@__DIR__,"create_MPC.jl"))
include(joinpath(@__DIR__,"iter_ref.jl"))
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
    ??::DELTA

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
        regularize!(KKT_reg, nx + ns,nz + ny, 1e-8) # TODO: update this ??
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
        ?? = DELTA(nx,ns,nz,ny)

        #initialization
        init = INIT(nx,nz,ny,idx,Q,G,A,idx_y .- ns)
        # nx,nz,ny,idx,Q,G,A,idx_y

        new(Q,q,A,b,G,h,x,s,z,y, KKT, KKT_reg, KKT_F, KKT_IR, rhs_a, rhs_c, p_a, p_c, idx, ??,CACHE(nx,nz,ny),init)
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
    # qp.??.x_a .= qp.p_a[qp.idx.x]
    qp.??.x_a .= view(qp.p_a,qp.idx.x)

    # qp.??.s_a .= qp.p_a[qp.idx.s]
    qp.??.s_a .= view(qp.p_a,qp.idx.s)

    # qp.??.z_a .= qp.p_a[qp.idx.z]
    qp.??.z_a .= view(qp.p_a,qp.idx.z)

    # qp.??.y_a .= qp.p_a[qp.idx.y]
    qp.??.y_a .= view(qp.p_a,qp.idx.y)
    return nothing
end
function index_sol_c!(qp::QP)
    qp.??.x_c .= view(qp.p_c,qp.idx.x)
    qp.??.s_c .= view(qp.p_c,qp.idx.s)
    qp.??.z_c .= view(qp.p_c,qp.idx.z)
    qp.??.y_c .= view(qp.p_c,qp.idx.y)
    return nothing
end
function linesearch(x,dx)
    # ?? = min(1.0, minimum([dx[i]<0 ? -x[i]/dx[i] : Inf for i = 1:length(x)]))
    ?? = 1.0
    for i = 1:length(x)
        if dx[i]<0
            ?? = min(??,-x[i]/dx[i])
        end
    end
    return ??
end
function centering_params(qp::QP)

    ?? = dot(qp.s,qp.z)/qp.idx.ns

    ?? = min(linesearch(qp.s,qp.??.s_a), linesearch(qp.z,qp.??.z_a))

    # ?? = (dot(qp.s + ??*qp.??.s_a, qp.z + ??*qp.??.z_a)/dot(qp.s,qp.z))^3
    @. qp.cache.z.c1 = qp.s + ??*qp.??.s_a
    @. qp.cache.z.c2 = qp.z + ??*qp.??.z_a
    ?? = (dot(qp.cache.z.c1,qp.cache.z.c2)/dot(qp.s,qp.z))^3

    return ??, ??
end

function rhs_kkt_c!(qp::QP, ??, ??)
    idx = qp.idx
    qp.rhs_c .= 0
    # qp.rhs_c[idx.s] = (??*?? .- (qp.??.s_a .* qp.??.z_a)) ./ qp.s
    @. qp.rhs_c[idx.s] = (??*?? - (qp.??.s_a * qp.??.z_a)) / qp.s
    return nothing
end
function combine_deltas!(qp::QP)
    @. qp.??.x = qp.??.x_a + qp.??.x_c
    @. qp.??.s = qp.??.s_a + qp.??.s_c
    @. qp.??.z = qp.??.z_a + qp.??.z_c
    @. qp.??.y = qp.??.y_a + qp.??.y_c
    return nothing
end
function update_vars!(qp::QP,??)
    @. qp.x += ??*qp.??.x
    @. qp.s += ??*qp.??.s
    @. qp.z += ??*qp.??.z
    @. qp.y += ??*qp.??.y
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
        ??, ?? = centering_params(qp)
        rhs_kkt_c!(qp, ??, ??)
        iterative_ref(qp.KKT_F, qp.KKT, qp.p_c, qp.rhs_c, qp.KKT_IR)
        index_sol_c!(qp)
        # @show norm(qp.KKT\qp.rhs_c - qp.p_c)
        # combine deltas
        combine_deltas!(qp)

        # last linesearch
        ?? = min(1,0.99*min(linesearch(qp.s,qp.??.s),linesearch(qp.z,qp.??.z)))

        update_vars!(qp,??)

        if logging(qp::QP,i,??)
            break
        end

    end
    return nothing
end


# # qp.rhs_a[idx.x] = -(qp.A'*qp.y + qp.G'*qp.z + qp.Q*qp.x + qp.q)
# mul!(c.x.c1,qp.A',qp.y)
# mul!(c.x.c2,qp.G',qp.z)
# mul!(c.x.c3,qp.Q,qp.x)
# @. qp.rhs_a[idx.x] = -c.x.c1 - c.x.c2 - c.x.c3 - qp.q
#
# # qp.rhs_a[idx.s] = -(qp.z)
# @. qp.rhs_a[idx.s] = -qp.z
#
# # qp.rhs_a[idx.z] = -(qp.G*qp.x + qp.s - qp.h)
# mul!(c.z.c1,qp.G,qp.x)
# @. qp.rhs_a[idx.z] = -c.z.c1 - qp.s + qp.h
#
# # qp.rhs_a[idx.y] = -(qp.A*qp.x - qp.b)
# mul!(c.y.c1,qp.A,qp.x)
# @. qp.rhs_a[idx.y] = -c.y.c1 + qp.b

function logging(qp::QP,iter,??)

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


    # eq_res = norm(qp.A*qp.x - qp.b)
    # ineq_res = norm(qp.G*qp.x + qp.s - qp.h)


    @printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter, J, gap, eq_res,
          ineq_res, ??)

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

    ?? = -Inf
    for i = 1:idx.ns
        if ?? < qp.z[i]
            ?? = qp.z[i]
        end
    end

    if ?? < 0
        for i = 1:idx.ns
            qp.s[i] = -qp.z[i]
        end
    else
        ?? += 1
        for i = 1:idx.ns
            qp.s[i] = -qp.z[i] + ??
        end
    end

    ?? = -Inf
    for i = 1:idx.ns
        if ?? < -qp.z[i]
            ?? = -qp.z[i]
        end
    end

    if ?? >= 0
        ?? += 1
        for i = 1:idx.ns
            qp.z[i] = qp.z[i] + ??
        end
    end


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
function regularize!(A,np,nd,??)
    for i = 1:np
        A[i,i] += ??
    end
    for i = (np + 1):(np + nd)
        A[i,i] -= ??
    end
    return nothing
end

@inline function update_kkt_factor!(qp::QP)
    QDLDL.update_A!(qp.KKT_F,qp.KKT_reg)
    return nothing
end
function update_kkt_reg!(qp::QP,??)
    idx = qp.idx
    # qp.KKT[idx.s, idx.s] = Diagonal(qp.z ./ qp.s)
    # qp.KKT[idx.s, idx.z] = I(idx.ns)
    for i = 1:idx.ns
        qp.KKT_reg[idx.s[i],idx.s[i]] = qp.KKT[idx.s[i],idx.s[i]] + ??
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
    N = 5
    Q,q,A,b,G,h = create_MPC(N)
    x1 = quadprog(Q,q,A,b,G,h)
    #
    # A = [A;G]
    # l = [b;-Inf*ones(length(h))]
    # u = [b;h]
    #
    # m = OSQP.Model()
    # OSQP.setup!(m; P = Q, q = q, A = A, l = l, u = u, eps_abs = 1e-8, eps_rel = 1e-8)
    #
    #
    # # @btime results = OSQP.solve!($m)
    # results = OSQP.solve!(m)
    # x = Variable(length(q))
    # problem = minimize(0.5*quadform(x,Matrix(Q)) + dot(q,x),[A*x == b, G*x <= h])
    # Convex.solve!(problem,Mosek.Optimizer)
    # @test norm(x.value - x1) <1e-3
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
