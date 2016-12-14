function jacobi(A::AbstractMatrix,b::AbstractVector;xinit=zero(b),tol=1e-10)
    m,n = size(A)
    @assert m == n
    x = xinit
    r = A*x-b
    D = Diagonal(diag(A))
    k=0
    while norm(r) > tol
        x = D\(b-(A*x-D*x))
        r = A*x-b
        k+=1
    end
    return x,k
end

function gaussseidel(A::AbstractMatrix,b::AbstractVector;xinit=zero(b),tol=1e-10)
    m,n = size(A)
    @assert m == n
    x = xinit
    r = A*x-b
    LPD = LowerTriangular(A)
    U = UpperTriangular(A)-Diagonal(diag(A))
    k=0
    while norm(r) > tol
        x = LPD\(b-U*x)
        r = A*x-b
        k+=1
    end
    return x,k
end

function pcg(A::AbstractMatrix,b::AbstractVector,P::AbstractMatrix;tol=1e-10)
	n = length(b)
	n1,n2 = size(A)
	n3,n4 = size(P)
	n == n1 == n2 == n3 == n4 || throw(DimensionMismatch(""))
	x = P\b
	r = b-A*x
	q = P\r
	p = q
	rq = r⋅q
    k=0
    while sqrt(abs(rq)) > tol
		Ap = A*p
		α = rq/(p⋅Ap)
		x = x+α*p
		r = r-α*Ap
		q = P\r
		rqnew = r⋅q
		p = q+rqnew/rq*p
		rq = rqnew
        k+=1
	end
    return x,k
end
