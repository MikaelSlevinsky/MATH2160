const math2160 = "/Users/Mikael/Documents/Courses/MATH2160"

using PyPlot

function computeHouseholder(A::Matrix,col::Int)
    u = A[:,col]
    u[1:col-1] = 0
    α = -copysign(sqrt(dot(u,u)), u[col])
    w = u
    w[col] = w[col]-α
    w
end

function applyHouseholder!(A::Matrix,w::Vector)
    # Psychologically H(w) = I - 2/⟨w,w⟩ ww^⊤,
    # but we implement this in-place on A.
    m,n = size(A)
    wTA = A'*w
    twoiww = 2/dot(w, w)
    scale!(twoiww, wTA)
    for j = 1:n
        @inbounds @simd for i = 1:m
            A[i,j] -= w[i]*wTA[j]
        end
    end
    A
end

function myQR(A::Matrix)
    R = copy(A)
    m,n = size(R)
    w = Array(Vector{Float64},n-1)
    for col = 1:n-1
        w[col] = computeHouseholder(R,col)
        applyHouseholder!(R,w[col])
    end
    Q = eye(max(m,n))
    for col = n-1:-1:1
        applyHouseholder!(Q,w[col])
    end
    Q,UpperTriangular(R)
end

include("IterativeLinearSolvers.jl")
include("SymTriPlusRankOne.jl")
