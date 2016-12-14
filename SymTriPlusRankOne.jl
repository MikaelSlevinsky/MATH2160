"""
    ..  SymTridiagonalPlusRankOne(T, u)

    Construct a real symmetric matrix `A` given by the addition of a symmetric tridiagonal matrix `T` and a symmetric rank-one outer product `uu^⊤`.

    A = T + uu^⊤

    The result is of type `SymTridiagonalPlusRankOne` and has an efficient matrix-vector product.
"""
type SymTridiagonalPlusRankOne <: AbstractMatrix{Float64}
    T::SymTridiagonal{Float64}
    u::Vector{Float64}
end

import Base: getindex, size, full, *, Diagonal, SymTridiagonal, LowerTriangular, UpperTriangular

getindex(A::SymTridiagonalPlusRankOne,i::Int,j::Int) = A.T[i,j]+A.u[i]*A.u[j]
getindex(A::SymTridiagonalPlusRankOne,i::Range,j::Int) = Float64[A[ii,j] for ii in i]
getindex(A::SymTridiagonalPlusRankOne,i::Int,j::Range) = Float64[A[i,jj] for jj in j].'
getindex(A::SymTridiagonalPlusRankOne,i::Range,j::Range) = Float64[A[ii,jj] for ii in i, jj in j]

size(A::SymTridiagonalPlusRankOne) = size(A.T)

full(A::SymTridiagonalPlusRankOne) = ((m,n)=size(A); A[1:m,1:n])

*(A::SymTridiagonalPlusRankOne,b::AbstractVector) = A.T*b + A.u.*(A.u.'*b)
*(A::SymTridiagonalPlusRankOne,b::AbstractMatrix) = A.T*b + A.u.*(A.u.'*b)

Diagonal(A::SymTridiagonalPlusRankOne) = Diagonal(diag(A))
SymTridiagonal(A::SymTridiagonalPlusRankOne) = SymTridiagonal(diag(A),[A[i-1,i] for i=2:size(A,1)])
LowerTriangular(A::SymTridiagonalPlusRankOne) = LowerTriangular(full(A))
UpperTriangular(A::SymTridiagonalPlusRankOne) = UpperTriangular(full(A))
