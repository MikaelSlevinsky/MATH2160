import Base: size, getindex, *, \

using LinearAlgebra

# We write a function which computes a vector `w` which parameterizes a
# Householder matrix H(w), so as to introduce zeros in the `col` column of `A`
# below the main diagonal. (See Lemma 2.3.8)

function computeHouseholder(A::Matrix, col::Integer)
    u = A[:, col]
    u[1:col-1] .= 0
    α = sign(u[col])*norm(u)
    w = u
    w[col] = α + w[col]
    w
end

# Next, we write a function which applies the Householder matrix to the matrix
# `A`. Notice that names of in-place functions are appended by `!`.
# (See Example 2.2.11)

function applyHouseholder!(A::Matrix, w::Vector)
    # Psychologically H(w) = I - 2/⟨w,w⟩ ww^⊤,
    # but we implement this in-place on A.
    m, n = size(A)
    wTA = transpose(w'A)
    twoiww = 2/w'w
    lmul!(twoiww, wTA)
    for j = 1:n
        @inbounds @simd for i = 1:m
            A[i,j] -= w[i]*wTA[j]
        end
    end
    A
end

"""
An unpivoted QR factorization for educational purposes.

# Examples
```jldoctest
julia> A = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4

julia> Q, R = myQR(A)
([-0.3162277660168382 -0.9486832980505137; -0.9486832980505138 0.3162277660168381], [-3.162277660168379 -4.427188724235732; 0.0 -0.632455532033676])

julia> norm(Q'Q - I) < norm(Q)*eps() # `Q` is orthogonal (subject to rounding errors)
true

julia> istriu(R) # `R` is upper triangular
true

julia> opnorm(Q) # See Lemma 2.2.8
1.0

julia> cond(Q) # See Theorem 2.2.9
1.0

julia> opnorm(A-Q*R) < sum(size(A))*eps()*opnorm(A) # The factorization has 2-normwise small residual relative to the 2-norm of `A`
true

julia> # For a computational exercise, try building your own QR factorization via Givens rotations.

```
"""
function myQR(A::AbstractMatrix)
    T = eltype(sqrt(one(eltype(A))))
    R = Matrix{T}(A)
    m, n = size(R)
    w = Vector{Vector{T}}(undef, n-1)
    for col = 1:n-1
        w[col] = computeHouseholder(R, col)
        applyHouseholder!(R, w[col])
    end
    Q = Matrix{T}(I, m, n)
    for col = n-1:-1:1
        applyHouseholder!(Q, w[col])
    end
    Q, UpperTriangular(R)
end

# Jacobi, Gauss-Seidel and pre-conditioned conjugate gradients.

function jacobi(A::AbstractMatrix, b::AbstractVector; xinit=zero(b), tol=1e-10)
    size(A, 1) == size(A, 2) == length(b) || throw(DimensionMismatch("Matrix has size $(size(A)) whereas right-hand side has length $(length(b))."))
    x = xinit
    r = A*x-b
    D = Diagonal(diag(A))
    k = 0
    while norm(r) > tol
        x = D\(b-(A*x-D*x))
        r = A*x-b
        k += 1
    end
    return x, k
end

function gaussseidel(A::AbstractMatrix, b::AbstractVector; xinit=zero(b), tol=1e-10)
    size(A, 1) == size(A, 2) == length(b) || throw(DimensionMismatch("Matrix has size $(size(A)) whereas right-hand side has length $(length(b))."))
    x = xinit
    r = A*x-b
    LPD = LowerTriangular(A)
    U = UpperTriangular(A)-Diagonal(diag(A))
    k = 0
    while norm(r) > tol
        x = LPD\(b-U*x)
        r = A*x-b
        k += 1
    end
    return x, k
end

function pcg(A::AbstractMatrix, b::AbstractVector; P = I, tol=1e-10)
	n1, n2 = size(A)
    if P isa UniformScaling
        n3, n4 = n1, n2
    else
        n3, n4 = size(P)
    end
	n1 == n2 == n3 == n4 == length(b) || throw(DimensionMismatch("Matrix has size $(size(A)) and preconditioner has size ($(n3), $(n4)) whereas right-hand side has length $(length(b))."))
	x = P\b
	r = b-A*x
	q = P\r
	p = q
	rq = r⋅q
    k = 0
    while sqrt(abs(rq)) > tol
		Ap = A*p
		α = rq/(p⋅Ap)
		x = x+α*p
		r = r-α*Ap
		q = P\r
		rqnew = r⋅q
		p = q+rqnew/rq*p
		rq = rqnew
        k += 1
	end
    return x, k
end

# A novel Julia struct for A = T + uuᵀ.

struct SymTridiagonalPlusRankOne{S <: Real} <: AbstractMatrix{S}
    T::SymTridiagonal{S}
    u::Vector{S}
    function SymTridiagonalPlusRankOne{S}(T::SymTridiagonal{S}, u::Vector{S}) where S <: Real
        if size(T, 1) ≠ length(u)
            throw(DimensionMismatch("Symmetric tridiagonal matrix has size $(size(T)) but symmetric rank-one outer-product has size ($(length(u)), $(length(u))). They should be the same."))
        end
        new{S}(T, u)
    end
end

"""
    SymTridiagonalPlusRankOne(T::SymTridiagonal{S}, u::Vector{S}) where S <: Real

Construct a real symmetric matrix `A` given by the addition of a symmetric tridiagonal matrix `T` and a symmetric rank-one outer product `uuᵀ`.

    A = T + uuᵀ

The result is of type `SymTridiagonalPlusRankOne` and has an efficient matrix-vector product.

# Examples
```jldoctest
julia> T = SymTridiagonal(2ones(5), ones(4))
5×5 SymTridiagonal{Float64,Array{Float64,1}}:
 2.0  1.0   ⋅    ⋅    ⋅
 1.0  2.0  1.0   ⋅    ⋅
  ⋅   1.0  2.0  1.0   ⋅
  ⋅    ⋅   1.0  2.0  1.0
  ⋅    ⋅    ⋅   1.0  2.0

julia> u = ones(5)
5-element Array{Float64,1}:
 1.0
 1.0
 1.0
 1.0
 1.0

julia> A = SymTridiagonalPlusRankOne(T, u)
5×5 SymTridiagonalPlusRankOne{Float64}:
 3.0  2.0  1.0  1.0  1.0
 2.0  3.0  2.0  1.0  1.0
 1.0  2.0  3.0  2.0  1.0
 1.0  1.0  2.0  3.0  2.0
 1.0  1.0  1.0  2.0  3.0

julia> A - T - u*u'
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> A*u
5-element Array{Float64,1}:
 8.0
 9.0
 9.0
 9.0
 8.0

julia> norm(A\\(A*u) - u)/norm(u)
4.550560269027491e-16

```
"""
SymTridiagonalPlusRankOne(T::SymTridiagonal{S}, u::Vector{S}) where S <: Real = SymTridiagonalPlusRankOne{S}(T, u)

size(A::SymTridiagonalPlusRankOne) = size(A.T)

getindex(A::SymTridiagonalPlusRankOne, i::Integer, j::Integer) = A.T[i,j] + A.u[i]*A.u[j]

*(A::SymTridiagonalPlusRankOne{S}, b::AbstractVector{U}) where {S <: Real, U} = A.T*b + A.u*(A.u'b)
*(A::SymTridiagonalPlusRankOne{S}, b::AbstractMatrix{U}) where {S <: Real, U} = A.T*b + A.u*(A.u'b)

# From the Sherman-Morrison formula (T+uuᵀ)⁻¹ = T⁻¹ - (T⁻¹uuᵀT⁻¹)/(1+uᵀT⁻¹u)
function \(A::SymTridiagonalPlusRankOne{S}, b::AbstractVector{U}) where {S <: Real, U}
    T = A.T
    u = A.u
    x = T\b
    Tiu = T\u
    x = x - inv(1+u'Tiu)*Tiu*(u'x)
    return x
end

function \(A::SymTridiagonalPlusRankOne{S}, b::AbstractMatrix{U}) where {S <: Real, U}
    T = A.T
    u = A.u
    x = T\b
    Tiu = T\u
    x = x - inv(1+u'Tiu)*Tiu*(u'x)
    return x
end
