### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ d33a4d92-fe0e-11ea-3cbf-c74e1daf820c
using LinearAlgebra

# ╔═╡ ad860ce4-fe0e-11ea-332f-dff64e52ed27
using Plots

# ╔═╡ 937d7990-fe0e-11ea-375c-698d6a5303b7
md"""Suppose we are interested in implementing a special linear algebra operation:

$$y = (A+A^\top)x + b.$$

In Julia, this is quite easy to do. First, we would check that $A$ is square (otherwise, the addition of $A$ to its transpose does not conform.). Then, we set the entries of $y$ by a nested pair of `for` loops given the fact that:

$$(Ax)_i = \sum_{j=1}^n A_{i,j} x_j\quad{\rm and}\quad (A^\top x)_i = \sum_{j=1}^nA_{j,i}x_j.$$"""

# ╔═╡ 9c4c979a-fe0e-11ea-3f56-c136cc0798bb
begin
	function my_special_problem!(y::Vector, A::Matrix, x::Vector, b::Vector)
		m, n = size(A)
		m ≠ n && throw(DimensionMismatch("Matrix A is not square."))
		length(y) == n == length(x) == length(b) || throw(DimensionMismatch("One of these vectors is behaving like an orangutan."))
		for i = 1:n
			y[i] = b[i]
			for j = 1:n
				y[i] += (A[i,j]+A[j,i])*x[j]
			end
		end
		return y
	end
	function my_special_problem2!(y::Vector, A::Matrix, x::Vector, b::Vector)
		m, n = size(A)
		m ≠ n && throw(DimensionMismatch("Matrix A is not square."))
		length(y) == n == length(x) == length(b) || throw(DimensionMismatch("One of these vectors is behaving like an orangutan."))
		for i = 1:n
			y[i] = b[i]
		end
		for i = 1:n
			for j = 1:n
				y[i] += A[j,i]*x[j]
				y[j] += A[j,i]*x[i]
			end
		end
		return y
	end
end

# ╔═╡ a07cec48-fe0e-11ea-1ac3-5592ff4f2859
y = zeros(15)

# ╔═╡ b55cbfe4-fe0e-11ea-1775-1791d0a29260
A = rand(15, 15)

# ╔═╡ b9f2e47a-fe0e-11ea-03ac-b1e1aaac6266
x = rand(15)

# ╔═╡ bd38a1e2-fe0e-11ea-21da-2f2464d7203f
b = rand(15)

# ╔═╡ bf5fac4a-fe0e-11ea-21b9-a5bd1ed3d0ed
my_special_problem!(y, A, x, b)

# ╔═╡ cc449ad8-fe0e-11ea-3c4d-23b183f6eb5d
norm(y - ((A+A')*x + b))

# ╔═╡ d0026c9a-fe0e-11ea-13db-8dd3a277d667
md"""When writing a special linear algebra problem such as the one above, it is important to keep track of the number of flops (floating-point operations), because this plays a basic role in prediciting the computational time. For any $n$, we can see by the nested `for` loops that after `y` is set to `b`, it takes $2n$ additions and $n$ multiplications for each $y_i$. In total, the function takes $3n^2$ flops. The second variant in fact costs $4n^2$ flops but is more friendly vis-à-vis the contiguous column-major storage of `A`.

When $n$ is small, the constants $3$ and $4$ may be helpful, but as $n$ gets large, what we really want to know is "if I double $n$, how much longer must I wait?" That question can be answered just with the exponent, $2$: a problem twice the size takes four times as long.

These arguments justify the new notation we are learning so we can summarize computational complexity more tersely and thus with more relevance: $3n^2 = \mathcal{O}(n^2)$ as $n\to\infty$."""

# ╔═╡ b026f0d0-fe0e-11ea-2a30-7362a62cb422
begin
	j = 1
	n = 2 .^ (1:12)
	t = zeros(length(n))
	r = zeros(length(n))
	for n in n
		y = zeros(n)
		A = rand(n, n)
		x = rand(n)
		b = rand(n)
		t[j] = @elapsed for i in 1:10 my_special_problem!(y, A, x, b) end
		t[j] /= 10
		r[j] = @elapsed for i in 1:10 my_special_problem2!(y, A, x, b) end
		r[j] /= 10
		global j += 1
	end
end

# ╔═╡ f5b23f60-fe0e-11ea-383a-17f83df994bb
begin
	scatter(n, t; xscale=:log10, yscale=:log10, label="my_special_problem!", legend=:topleft)
	scatter!(n, r; xscale=:log10, yscale=:log10, label="my_special_problem2!", legend=:topleft)
	plot!(n, 3e-9n.^2, label="\$\\mathcal{O}(n^2)\$")
	xlabel!("n")
	ylabel!("Execution Time (s)")
end

# ╔═╡ 1f40a46e-fe77-11ea-2717-ddf799cd567f
md"""In the early stages, the first implementation is a bit faster than the second. But as soon as memory complexity begins to play a role, respecting the contiguous ordering of `A` improves the timings by about a factor of four."""

# ╔═╡ 522c478c-fe77-11ea-0a65-7549475bd133
t[end]/r[end]

# ╔═╡ Cell order:
# ╟─937d7990-fe0e-11ea-375c-698d6a5303b7
# ╠═9c4c979a-fe0e-11ea-3f56-c136cc0798bb
# ╠═a07cec48-fe0e-11ea-1ac3-5592ff4f2859
# ╠═b55cbfe4-fe0e-11ea-1775-1791d0a29260
# ╠═b9f2e47a-fe0e-11ea-03ac-b1e1aaac6266
# ╠═bd38a1e2-fe0e-11ea-21da-2f2464d7203f
# ╠═bf5fac4a-fe0e-11ea-21b9-a5bd1ed3d0ed
# ╠═d33a4d92-fe0e-11ea-3cbf-c74e1daf820c
# ╠═cc449ad8-fe0e-11ea-3c4d-23b183f6eb5d
# ╟─d0026c9a-fe0e-11ea-13db-8dd3a277d667
# ╠═b026f0d0-fe0e-11ea-2a30-7362a62cb422
# ╠═ad860ce4-fe0e-11ea-332f-dff64e52ed27
# ╠═f5b23f60-fe0e-11ea-383a-17f83df994bb
# ╟─1f40a46e-fe77-11ea-2717-ddf799cd567f
# ╠═522c478c-fe77-11ea-0a65-7549475bd133
