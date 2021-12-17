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
