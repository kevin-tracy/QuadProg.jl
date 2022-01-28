using LinearAlgebra
using ForwardDiff
using Printf

include(joinpath(@__DIR__,"alloc_free/create_mpc.jl"))

const FD = ForwardDiff

struct IDX
    # contains the variable indexing stuff for β = [x;s;z;y]
    nx::Int64
    ns::Int64
    nz::Int64
    ny::Int64

    x::UnitRange{Int64}
    s::UnitRange{Int64}
    z::UnitRange{Int64}
    y::UnitRange{Int64}

    function IDX(q,h,b)
        # create indexing stuff from problem size
        nx = length(q)
        ns = length(h)
        nz = ns
        ny = length(b)
        new(nx,ns,nz,ny,
            1:nx,
            (nx + 1) : (nx + ns),
            (nx + ns + 1) : (nx + ns + nz),
            (nx + ns + nz + 1) : (nx + ns + nz + ny))
    end
end
struct QP
    # problem data
    Q::Matrix{Float64}
    q::Array{Float64,1}
    A::Matrix{Float64}
    b::Array{Float64,1}
    G::Matrix{Float64}
    h::Array{Float64,1}
end

function unpackβ(β,idx)
    # unpack the β = [x;s;z;y]
    x = β[idx.x]
    s = β[idx.s]
    z = β[idx.z]
    y = β[idx.y]
    return x,s,z,y
end

function kkt_conditions(qp,idx,β)
    x,s,z,y = unpackβ(β,idx)
    return [(qp.A'*y + qp.G'*z + qp.Q*x + qp.q); # stationarity
            (s .* z);                            # complimentarity
            (qp.G*x + s - qp.h);                 # primal feasibility
            (qp.A*x - qp.b)]                     # primal feasibility
end

function linesearch(x,dx)
    # this returns the max α ∈ [0,1] st (x + Δx > 0)
    α = min(1.0, minimum([dx[i] < 0 ? -x[i]/dx[i] : Inf for i = 1:length(x)]))
    return α
end
function centering_params(s,z,s_a,z_a)
    # surrogate duality gap
    μ = dot(s,z)/length(s)

    # find α such that s and z remain positive
    α = min(linesearch(s,s_a), linesearch(z,z_a))

    # Mehrotra
    σ = (dot(s + α*s_a, z + α*z_a)/dot(s,z))^3
    return σ, μ
end
function lazy_pdip(Q,q,A,b,G,h)

    # build QP
    qp = QP(Q,q,A,b,G,h)

    # create indices
    idx = IDX(q,h,b)

    # initial guess (s,z>0)
    x = zeros(idx.nx)
    s = ones(idx.ns)
    z = ones(idx.nz)
    y = zeros(idx.ny)
    β = [x;s;z;y]

    @printf "iter  gap        |Ax-b|     |Gx+s-h|   α\n"
    @printf "---------------------------------------------\n"

    # main loop
    for i = 1:30

        # evaluate kkt conditions
        r = kkt_conditions(qp,idx,β)

        # check convergence
        if norm(r,Inf)<1e-10
            break
        end

        # get jacobian of kkt system
        K = FD.jacobian(_β -> kkt_conditions(qp,idx,_β),β)

        # solve for affine Newton step
        Δβ_a = -K\r

        # centering + correcting
        σ, μ = centering_params(β[idx.s],β[idx.z],Δβ_a[idx.s],Δβ_a[idx.z])

        # solve for centering + correcting step
        r_c = [zeros(idx.nx);
               σ*μ .- (Δβ_a[idx.s] .* Δβ_a[idx.z]);
               zeros(idx.nz + idx.ny)]
        Δβ_c = K\r_c

        # combine steps
        Δβ = Δβ_a + Δβ_c

        # final line search
        α = min(1,0.99*min(linesearch(β[idx.s],Δβ[idx.s]),linesearch(β[idx.z],Δβ[idx.z])))

        # update β
        β += α*Δβ

        @printf("%3d  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          i, sum(r[idx.s])/idx.ns, norm(r[idx.y]),
          norm(r[idx.z]), α)

    end

end


function tt()

    Q,q,A,b,G,h = create_MPC(10)

    lazy_pdip(Q,q,A,b,G,h)

    return nothing

end

tt()
