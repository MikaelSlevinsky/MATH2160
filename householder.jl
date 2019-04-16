# This file uses three functions to write a QR factorization of a matrix
# from scratch. First, we load a standard library, `LinearAlgebra`, which
# ships with Julia v0.7 and v1.0.

using LinearAlgebra

# Then, we write a function which computes a vector `w` which parameterizes a
# Householder matrix H(w), so as to introduce zeros in the `col` column of `A`
# below the main diagonal. (See Lemma 2.3.8)

function computeHouseholder(A::Matrix, col::Int)
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
    lmul!(2/dot(w, w), wTA)
    for j = 1:n
        for i = 1:m
            A[i, j] -= w[i]*wTA[j]
        end
    end
    A
end

# Putting it all together, we have a simple function which loops through columns
# 1:n-1, creates the Householder reflectors that introduce the appropriate zeros
# in `A` (copied into `R`), and accumulates the Householder reflector in an
# orthogonal matrix `Q`. (See Theorems 2.3.6 & 2.2.12)

function myQR(A::Matrix)
    T = eltype(sqrt(one(eltype(A))))
    m, n = size(A)
    R = Matrix{T}(undef, m, n)
    for j = 1:n, i = 1:m
        R[i, j] = A[i, j]
    end
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

# Typical usage. Take any matrix `A`:

A = [1 2; 3 4]

# Calculate orthogonal and upper triangular factors:

Q, R = myQR(A)

# Check that `Q` is indeed orthogonal (subject to rounding errors on the order
# of machine precision):

opnorm(Q'Q - I) < sum(size(Q))*eps()

# See Lemma 2.2.8

opnorm(Q)

# See Theorem 2.2.9

cond(Q)

# Check that `R` is upper triangular:

istriu(R)

# Check that the factorization has 2-normwise small residual, relative to the
# 2-norm of the matrix itself.

opnorm(A-Q*R)/opnorm(A) < sum(size(A))*eps()


# For a computational exercise, try building your own QR factorization via Givens rotations.
