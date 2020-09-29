### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 582c51b8-01f6-11eb-0051-c9d04a4baf76
using LinearAlgebra, Plots

# ╔═╡ bb2e7882-01fa-11eb-3f9f-bf05ba047c4e
md"""Let's make and plot an $n$-sided random polygon with $2$-normalized vertices."""

# ╔═╡ 66a17c32-01f6-11eb-3a77-4bc6bbd981e9
begin
	n = 10
	x = normalize!(rand(n))
	y = normalize!(rand(n))
end

# ╔═╡ be026478-01f6-11eb-04a4-bb36e1e45ecf
begin
	scatter(x, y; legend=false)
	plot!([x; x[1]], [y; y[1]])
end

# ╔═╡ dd83c1da-01fa-11eb-155c-5b7b7d7e55b4
md"""What happens when we consider the sequence of polygons whose vertices are the midpoints of previous polygon's vertices?"""

# ╔═╡ dc3a6d6e-01f6-11eb-1cdb-b300919a2ac4
begin
	x1 = x[1]
	y1 = y[1]
	for i in 1:length(x)-1
		x[i] = (x[i]+x[i+1])/2
		y[i] = (y[i]+y[i+1])/2
	end
	x[end] = (x[end]+x1)/2
	y[end] = (y[end]+y1)/2
	scatter!(x, y; legend=false)
	plot!([x; x[1]], [y; y[1]])
end

# ╔═╡ c3cfdd5e-01f6-11eb-0c11-7f10d33c40f1
md"""**Observation 1:** they averaged polygons converge to a point, the centroid:

$$\bar{x} = \frac{1}{n}\sum_{i=1}^nx_i,\quad{\rm and}\quad \bar{y} = \frac{1}{n}\sum_{i=1}^n x_i.$$
"""

# ╔═╡ 8199ca4e-01fb-11eb-079a-cfbddc175166
sum(x)/n, sum(y)/n

# ╔═╡ ce83c5d0-01fb-11eb-1840-8519560f174b
md"""If we work with centroid-$0$ polygons and renormalize at every iteration, we will no longer converge to a point, but perhaps we can observe other phenomena."""

# ╔═╡ 01a8f9b2-01fc-11eb-1a3a-e3c5f7212527
begin
	m = 25
	u = normalize!(rand(m))
	v = normalize!(rand(m))
	u .-= sum(u)/m
	v .-= sum(v)/m
end

# ╔═╡ 0f66546e-01fc-11eb-0ed2-dd693dbe00f0
begin
	scatter(u, v; legend=false)
	plot!([u; u[1]], [v; v[1]])
end

# ╔═╡ 3058fd3e-01fc-11eb-03c2-195f7318b802
begin
	u1 = u[1]
	v1 = v[1]
	for i in 1:length(u)-1
		u[i] = (u[i]+u[i+1])/2
		v[i] = (v[i]+v[i+1])/2
	end
	u[end] = (u[end]+u1)/2
	v[end] = (v[end]+v1)/2
	normalize!(u)
	normalize!(v)
	scatter(u, v; legend=false)
	plot!([u; u[1]], [v; v[1]])
end

# ╔═╡ 824b65dc-01fc-11eb-1ee3-15c67bfc969b
md"""**Observation 2:** The vertices of the polygons converge to the boundary of an ellipse at a $45^\circ$ angle."""

# ╔═╡ ace4bcc6-01fc-11eb-077b-27c79a6eaafc
md"""**Observation 3:** After converging (to plotting accuracy), polygons of even iterates (and odd iterates) appear almost the same."""

# ╔═╡ 121f1744-01fd-11eb-1727-1fbcd41649b1
md"""The averaging of vertices can be described in terms of a matrix-vector product:

$$\begin{bmatrix} x_1^{(k+1)}\\ x_2^{(k+1)}\\ \vdots\\ x_{n-1}^{(k+1)}\\ x_n^{(k+1)} \end{bmatrix} = \frac{1}{2}\begin{bmatrix} 1 & 1\\ & 1 & 1\\ & & \ddots & \ddots\\ & & & \ddots & 1\\ 1 & & & & 1\end{bmatrix} \begin{bmatrix} x_1^{(k)}\\ x_2^{(k)}\\ \vdots\\ x_{n-1}^{(k)}\\ x_n^{(k)} \end{bmatrix},$$
with a similar equation holding for $y$.

This matrix can be described in terms of the $n\times n$ identity, $I_n$, and the circular shift:

$$S_n = \begin{bmatrix} 0 & 1\\ & 0 & 1\\ & & \ddots & \ddots\\ & & & \ddots & 1\\ 1 & & & & 0\end{bmatrix},$$

so that $x^{(k+1)} = \frac{1}{2}(I_n+S_n)x^{(k)}.$

In fact, not a single entry needs to be stored in order to apply this matrix to a vector. In Julia, the *array interface* consists of at least two functions: `size` and `getindex`."""

# ╔═╡ 5023a522-01fe-11eb-1b42-195e77262398
begin
	struct HalfIdentityPlusCircularShift{T} <: AbstractMatrix{T}
		n::Int
	end
	Base.size(A::HalfIdentityPlusCircularShift{T}) where T = (A.n, A.n)
	function Base.getindex(A::HalfIdentityPlusCircularShift{T}, i, j) where T
		n = A.n
		if (i == j || i+1 == j) && 1 ≤ i ≤ n && 1 ≤ j ≤ n
			return T(0.5)
		elseif i == n && j == 1
			return T(0.5)
		else
			return T(0)
		end
	end
end

# ╔═╡ 56ba75a4-01ff-11eb-3ec9-eda030f81717
HalfIdentityPlusCircularShift{Float64}(n)

# ╔═╡ ce8463b0-01ff-11eb-1774-515e70e26a1b
HalfIdentityPlusCircularShift{Rational{Int}}(n)

# ╔═╡ 93d86806-0263-11eb-0812-2bbef4d4be60
md"""By subtyping our struct as an `AbstractMatrix` and by implementing `size` and `getindex`, we get certain extras for free."""

# ╔═╡ 8e3f8da2-0263-11eb-01c8-afc9eadf49c2
HalfIdentityPlusCircularShift{Float64}(n)*x

# ╔═╡ f4f44470-01ff-11eb-0216-d1bc663b60da
md"""For more information on the random polygons, please see

	1. A. N. Elmachtoub and C. F. Van Loan, [From Random Polygon to Ellipse: An Eigenanalysis](https://epubs.siam.org/doi/pdf/10.1137/090746707), *SIAM Rev.*, **52**:151--170, 2010."""

# ╔═╡ Cell order:
# ╠═582c51b8-01f6-11eb-0051-c9d04a4baf76
# ╟─bb2e7882-01fa-11eb-3f9f-bf05ba047c4e
# ╠═66a17c32-01f6-11eb-3a77-4bc6bbd981e9
# ╠═be026478-01f6-11eb-04a4-bb36e1e45ecf
# ╟─dd83c1da-01fa-11eb-155c-5b7b7d7e55b4
# ╠═dc3a6d6e-01f6-11eb-1cdb-b300919a2ac4
# ╟─c3cfdd5e-01f6-11eb-0c11-7f10d33c40f1
# ╠═8199ca4e-01fb-11eb-079a-cfbddc175166
# ╟─ce83c5d0-01fb-11eb-1840-8519560f174b
# ╠═01a8f9b2-01fc-11eb-1a3a-e3c5f7212527
# ╠═0f66546e-01fc-11eb-0ed2-dd693dbe00f0
# ╠═3058fd3e-01fc-11eb-03c2-195f7318b802
# ╟─824b65dc-01fc-11eb-1ee3-15c67bfc969b
# ╟─ace4bcc6-01fc-11eb-077b-27c79a6eaafc
# ╟─121f1744-01fd-11eb-1727-1fbcd41649b1
# ╠═5023a522-01fe-11eb-1b42-195e77262398
# ╠═56ba75a4-01ff-11eb-3ec9-eda030f81717
# ╠═ce8463b0-01ff-11eb-1774-515e70e26a1b
# ╟─93d86806-0263-11eb-0812-2bbef4d4be60
# ╠═8e3f8da2-0263-11eb-01c8-afc9eadf49c2
# ╟─f4f44470-01ff-11eb-0216-d1bc663b60da
