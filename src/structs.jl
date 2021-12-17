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
