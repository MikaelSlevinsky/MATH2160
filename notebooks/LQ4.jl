### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 9cccb8a6-0351-11eb-2afe-63c84e8519e1
using LinearAlgebra

# ╔═╡ eb5da5fe-034f-11eb-15db-23a210da52db
md"""In the [RandomPolygons](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/RandomPolygons.jl.html) notebook, we created a random polygon and new polygons by averaging the old polygon's vertices. The averaging of vertices can be described in terms of a matrix-vector product:

$$\begin{bmatrix} x_1^{(k+1)}\\ x_2^{(k+1)}\\ \vdots\\ x_{n-1}^{(k+1)}\\ x_n^{(k+1)} \end{bmatrix} = \frac{1}{2}\begin{bmatrix} 1 & 1\\ & 1 & 1\\ & & \ddots & \ddots\\ & & & \ddots & 1\\ 1 & & & & 1\end{bmatrix} \begin{bmatrix} x_1^{(k)}\\ x_2^{(k)}\\ \vdots\\ x_{n-1}^{(k)}\\ x_n^{(k)} \end{bmatrix},$$
with a similar equation holding for $y$.

This matrix can be described in terms of the $n\times n$ identity, $I_n$, and the circular shift:

$$S_n = \begin{bmatrix} 0 & 1\\ & 0 & 1\\ & & \ddots & \ddots\\ & & & \ddots & 1\\ 1 & & & & 0\end{bmatrix},$$

so that $x^{(k+1)} = \frac{1}{2}(I_n+S_n)x^{(k)}.$

In fact, not a single entry needs to be stored in order to apply this matrix to a vector. We implemented the most basic aspects of Julia's array interface by coding methods for `size` and `getindex`:"""

# ╔═╡ d751747c-0350-11eb-3a0e-75cf0161be8d
begin
	struct HalfIdentityPlusCircularShift{T} <: AbstractMatrix{T}
		n::Int
	end
	Base.size(A::HalfIdentityPlusCircularShift{T}) where T = (A.n, A.n)
	function Base.getindex(A::HalfIdentityPlusCircularShift{T}, i::Int, j::Int) where T
		n = A.n
		if (i == j || i+1 == j) && 1 ≤ i ≤ n && 1 ≤ j ≤ n
			return T(0.5)
		elseif i == n && j == 1
			return T(0.5)
		else
			return T(0)
		end
	end
	struct HalfIdentityPlusCircularShiftFactors{T} <: AbstractMatrix{T}
    	n::Int
	end
	Base.size(A::HalfIdentityPlusCircularShiftFactors{T}) where T = (A.n, A.n)
	function Base.getindex(A::HalfIdentityPlusCircularShiftFactors{T}, i::Int, j::Int) where T
		# Insert correct indexing here, Read to the end before coming back.
		return T(0)
	end
end

# ╔═╡ 0d57a88e-0351-11eb-297d-45397dcd4938
A = HalfIdentityPlusCircularShift{Float64}(11)

# ╔═╡ f5ac075c-0350-11eb-30bf-4beb7284a156
md"""The array interface gives us a lot of functionality for *free*, even if it is not the most efficient implementation. In the next block, check that you can compute:
 - the square of the matrix, $A^2$;
 - the product $A^\top A$;
 - the trace, ${\rm tr}(A)$; and,
 - the eigenvalues of $A$ (by `eigvals`).
Do not change the names of the variables.
"""

# ╔═╡ c489d446-038d-11eb-041b-076c0c62bd3e
A2 = A # Wrong! Fix me.

# ╔═╡ c4484bde-038d-11eb-291b-cff93ef45ab7
AtA = A # Wrong! Fix me.

# ╔═╡ c3c6fdcc-038d-11eb-1292-8b0b6626261c
trA = 0.0 # Wrong! Fix me.

# ╔═╡ c329ef5a-038d-11eb-27f3-e7db3593d808
λ = ones(10) # Wrong! Fix me.

# ╔═╡ b5409dd0-035b-11eb-2762-edce74db7f48
if !(@isdefined A2) || !(@isdefined AtA) || !(@isdefined trA) || !(@isdefined λ)
	md"""Do not change the name of the variables - write your answer as
```
A2 = "..."
AtA = "..."
trA = "..."
λ = "..."
```
	"""
end

# ╔═╡ 67383d4a-0352-11eb-3199-0deba4ad424b
let A = Matrix{Float64}(HalfIdentityPlusCircularShift{Float64}(11))
	if A2 ≈ A^2 && AtA ≈ A'A && trA ≈ tr(A) && λ ≈ eigvals(A)
		md"""
!!! correct
    Well done!
		"""
	else
		md"""
!!! warning "Incorrect"
    Keep working on it!
		"""
	end
end

# ╔═╡ 4ce5233a-0353-11eb-0cc2-4df1cda7d738
md"""A matrix-vector product with $A = \frac{1}{2}(I_n+S_n)$ costs only $\mathcal{O}(n)$ floating-point operations (flops) rather than the direct $\mathcal{O}(n^2)$ flops. Let's define a new method for `LinearAlgebra.mul!` to reset $y\leftarrow Ax$. When we call `A*x`, Julia will magically know to use this new method."""

# ╔═╡ f2e5ea12-0350-11eb-15bd-599472b0d1f4
begin
	function LinearAlgebra.mul!(y::Vector{T}, A::HalfIdentityPlusCircularShift{T}, x::Vector{T}) where T
		n = size(A, 1)
		if n != length(x)
			throw(DimensionMismatch("matrix A has dimensions $(size(A)), vector x has length $(length(x))"))
		end
		if n != length(y)
			throw(DimensionMismatch("result y has length $(length(y)), needs length $n."))
		end
		# Insert floating-point arithmetic here
		return y
	end
	x = collect(1.0:11)
	y = A*x
end

# ╔═╡ 8de52df0-0353-11eb-1bf7-89f6b6a6b60f
let A = Matrix{Float64}(HalfIdentityPlusCircularShift{Float64}(11))
	if y ≈ A*x
		md"""
!!! correct
    Well done!
		"""
	else
		md"""
!!! warning "Incorrect"
    Keep working on it!
		"""
	end
end

# ╔═╡ 14f8e460-0354-11eb-2af6-ef45d2501911
md"""We can also improve upon the solution of linear systems involving $A$, given the fact that the unpivoted $LU$ decomposition is particularly easy to write down based on:

$$\frac{1}{2}\begin{bmatrix} 1 & 1\\ & 1 & 1\\ & & \ddots & \ddots\\ & & & \ddots & 1\\ 1 & & & & 1\end{bmatrix} = \frac{1}{2}\begin{bmatrix} 1\\ & 1\\ & & 1\\ & & & \ddots\\ 1 & -1 & 1 & \cdots & 1\end{bmatrix} \begin{bmatrix} 1 & 1\\ & 1 & 1\\ & & \ddots & \ddots\\ & & & 1 & 1\\ & & & & 1-(-1)^n\end{bmatrix}.$$

Go back to the top and add to the second cell the missing code in the `getindex` method for the type `HalfIdentityPlusCircularShiftFactors`. This will return the entries of $L$ strictly below the main diagonal and $U$ on the main diagonal and above.
"""

# ╔═╡ 00cc07cc-0356-11eb-22ba-a1b62ab96161
begin
	function LinearAlgebra.lu(A::HalfIdentityPlusCircularShift{T}, pivot::Val{false}=Val(false); check::Bool=true) where T
		n = size(A, 1)
		check && iseven(n) && throw(SingularException(n))
		F = HalfIdentityPlusCircularShiftFactors{T}(n)
		LU{T, typeof(F)}(F, collect(1:n), isodd(n) ? 0 : 1)
	end
	function Base.getproperty(F::LU{T,HalfIdentityPlusCircularShiftFactors{T}}, d::Symbol) where T
		m, n = size(F)
		if d === :L
			L = tril!(getfield(F, :factors)[1:m, 1:min(m,n)])
			for i = 1:min(m,n); L[i,i] = one(T); end
			return L
		elseif d === :U
			return triu!(getfield(F, :factors)[1:min(m,n), 1:n])
		elseif d === :p
			return ipiv2perm(getfield(F, :ipiv), m)
		elseif d === :P
			return Matrix{T}(I, m, m)[:,invperm(F.p)]
		else
			getfield(F, d)
		end
	end
end

# ╔═╡ 5e0a7bfe-035f-11eb-0dbe-7724b6a2d914
F = lu(A)

# ╔═╡ 7b0726e4-035f-11eb-18ef-55f2b8d98bda
let A = Matrix{Float64}(A)
	G = lu(A, Val(false); check = false)
	if F.L ≈ G.L && F.U ≈ G.U
		md"""
!!! correct
    Well done!
		"""
	else
		md"""
!!! warning "Incorrect"
    Keep working on it!
		"""
	end
end

# ╔═╡ 7fefda8a-0359-11eb-00f0-ff6fdde6b6a8
md"""In the cell below, add the floating-point arithmetic that is needed to solve linear systems with the $LU$ decomposition of $A$."""

# ╔═╡ 396cb636-0358-11eb-1a07-51f60ff41760
begin
	function LinearAlgebra.ldiv!(F::LU{T,HalfIdentityPlusCircularShiftFactors{T}}, x::Vector{T}) where T
		n = size(F, 1)
		if n != length(x)
			throw(DimensionMismatch("matrix A has dimensions $(size(A)), vector x has length $(length(x))"))
		end
		# Insert floating-point arithmetic here
		return x
	end
	b = 1 ./ (1:11)
	z = A\b
end

# ╔═╡ 6c987128-0358-11eb-17ee-1788275a60e2
let A = Matrix{Float64}(A), b = 1 ./ (1:11)
	if z ≈ A\b
		md"""
!!! correct
    Well done!
		"""
	else
		md"""
!!! warning "Incorrect"
    Keep working on it!
		"""
	end
end

# ╔═╡ b7b9c280-0359-11eb-0297-678fa5497c41
md"""Describe, in your own words, why $A$ is singular for even $n$ and invertible when $n$ is odd."""

# ╔═╡ c807c53a-0359-11eb-27a2-f71b3fc4b2d1
paragraph = md"""Your answer here"""

# ╔═╡ fa5898fc-0359-11eb-024e-75b5b60d5306
if !(@isdefined paragraph)
	md"""Do not change the name of the variable - write your answer as `paragraph = "..."`"""
end

# ╔═╡ 6819a6e4-035d-11eb-120c-75a8ea369518
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ 628796ac-035d-11eb-0cc6-ad4e72554d81
hint(md"To reset $x\leftarrow A^{-1}x$ with an $LU$ factorization of $A$, first reset $x \leftarrow L^{-1}x$ by forward substitution, then reset $x \leftarrow U^{-1} x$ by backward substitution. Be sure to use only the necessary arithmetic.")

# ╔═╡ Cell order:
# ╟─eb5da5fe-034f-11eb-15db-23a210da52db
# ╠═d751747c-0350-11eb-3a0e-75cf0161be8d
# ╠═0d57a88e-0351-11eb-297d-45397dcd4938
# ╠═9cccb8a6-0351-11eb-2afe-63c84e8519e1
# ╟─f5ac075c-0350-11eb-30bf-4beb7284a156
# ╠═c489d446-038d-11eb-041b-076c0c62bd3e
# ╠═c4484bde-038d-11eb-291b-cff93ef45ab7
# ╠═c3c6fdcc-038d-11eb-1292-8b0b6626261c
# ╠═c329ef5a-038d-11eb-27f3-e7db3593d808
# ╟─b5409dd0-035b-11eb-2762-edce74db7f48
# ╟─67383d4a-0352-11eb-3199-0deba4ad424b
# ╟─4ce5233a-0353-11eb-0cc2-4df1cda7d738
# ╠═f2e5ea12-0350-11eb-15bd-599472b0d1f4
# ╟─8de52df0-0353-11eb-1bf7-89f6b6a6b60f
# ╟─14f8e460-0354-11eb-2af6-ef45d2501911
# ╠═00cc07cc-0356-11eb-22ba-a1b62ab96161
# ╠═5e0a7bfe-035f-11eb-0dbe-7724b6a2d914
# ╟─7b0726e4-035f-11eb-18ef-55f2b8d98bda
# ╟─7fefda8a-0359-11eb-00f0-ff6fdde6b6a8
# ╟─628796ac-035d-11eb-0cc6-ad4e72554d81
# ╠═396cb636-0358-11eb-1a07-51f60ff41760
# ╟─6c987128-0358-11eb-17ee-1788275a60e2
# ╟─b7b9c280-0359-11eb-0297-678fa5497c41
# ╠═c807c53a-0359-11eb-27a2-f71b3fc4b2d1
# ╟─fa5898fc-0359-11eb-024e-75b5b60d5306
# ╟─6819a6e4-035d-11eb-120c-75a8ea369518
