### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0fbb2042-f695-11ea-302e-f5c537a1e13e
using LinearAlgebra, Plots

# ╔═╡ 209cac78-f695-11ea-2f5a-a7730093e80c
md"""It's usually quite straightforward how to make sense of the vector $p$-norms:

$$\|x\|_1 = \sum_{i=1}^n|x_i|,\qquad \|x\|_2 = \sqrt{\sum_{i=1}^n|x_i|^2},\qquad \|x\|_\infty = \max_{1\le i\le n} |x_i|.$$

What is perhaps not so obvious are induced matrix norms:

$$\|A\|_p = \sup_{0\ne x\in\mathbb{C}^n} \dfrac{\|Ax\|_p}{\|x\|_p} = \max_{\|x\|_p\le1} \|Ax\|_p = \max_{\|x\|_p = 1} \|Ax\|_p.$$

To visualize these, we will restrict our attention to $m=n=2$; that is, vectors live in the plane. For $p=1$ and $p=\infty$, we get lozenges, and for $p=2$, we get ellipses.
"""

# ╔═╡ 0622a5ea-f696-11ea-3a9f-07f16326b72c
@bind θ html"""<input type="range" min=0 max=2 step=0.01 value=0>"""

# ╔═╡ 6c96e420-f697-11ea-32d5-7138a4a8df22
@bind p html"""<input type="range" min=1 max=10 step=0.1 value=1>"""

# ╔═╡ 05e09f60-f696-11ea-3084-27e6463d8ccd
θ, p

# ╔═╡ 05cb8ec2-f696-11ea-32a9-5339bece53b8
A = [1 2; 0 2.0]

# ╔═╡ 05aa91c2-f696-11ea-0cf7-1f20f038c9c5
begin
	t = 0:0.01:2
	v = [normalize!([cospi(t), sinpi(t)], p) for t in t]
	u = [A*v for v in v]
	plot(map(u->u[1], u), map(u->u[2], u); label="Map of the unit circle in the $p-norm.", legend = :topleft)
	x = normalize!([cospi(θ), sinpi(θ)], p)
	y = A*x
	plot!([0,y[1]], [0,y[2]]; label="Vector with angle $θ*π.")
	xlims!(-3,3)
	ylims!(-3,3)
end

# ╔═╡ b5a6b5b2-f69e-11ea-2b7b-ab3690d0a779
md"""For matrices, the induced $p$-norms are accessed in Julia via `opnorm`, short for operator norm (because a matrix is an archetypal linear operator)."""

# ╔═╡ d7e9bc64-f69e-11ea-0da3-2b4bf670defb
opnorm(A, 1), opnorm(A, 2), opnorm(A, Inf)

# ╔═╡ 5027ca24-f698-11ea-36c0-01e10e6a3677
md"""Similarly, it is important to grasp a strong intuition on $p$-norms of functions. This is most easily visualized when looking at an approximation error. Suppose we wish to approximation $f(x) = e^x$ on $[-1,1]$. Then a degree-$2$ polynomial interpolant through the points $\{-1,0,1\}$ has the form:

$$p_2(x) = e^{-1}\frac{x(x-1)}{2} - (x-1)(x+1) + e^1\frac{x(x+1)}{2}.$$

Let's plot the integrand in the $p$-norm:

$$\|f-p_2\|_p = \left(\int_{-1}^1 |f(x)-p_2(x)|^p{\rm\,d}x\right)^{\frac{1}{p}}.$$"""

# ╔═╡ 5042ab46-f698-11ea-3847-4b6867cdaa18
begin
	s = range(-1, stop = 1, length=1001)
	f = x -> exp(x)
	p2 = x -> exp(-1)*x*(x-1)/2 - (x-1)*(x+1) + exp(1)*x*(x+1)/2
	l = @layout([a b])
	plot(
		plot(s, [f.(s), p2.(s)]; label=permutedims(["\$f(x)\$", "\$p_2(x)\$"]), legend=:topleft),
		plot(s, abs.(f.(s) .- p2.(s)).^p; label="\$|f(x)-p_2(x)|^{$(p)}\$", legend=:topleft)
	)
end

# ╔═╡ Cell order:
# ╠═0fbb2042-f695-11ea-302e-f5c537a1e13e
# ╟─209cac78-f695-11ea-2f5a-a7730093e80c
# ╠═0622a5ea-f696-11ea-3a9f-07f16326b72c
# ╠═6c96e420-f697-11ea-32d5-7138a4a8df22
# ╠═05e09f60-f696-11ea-3084-27e6463d8ccd
# ╠═05cb8ec2-f696-11ea-32a9-5339bece53b8
# ╠═05aa91c2-f696-11ea-0cf7-1f20f038c9c5
# ╟─b5a6b5b2-f69e-11ea-2b7b-ab3690d0a779
# ╠═d7e9bc64-f69e-11ea-0da3-2b4bf670defb
# ╟─5027ca24-f698-11ea-36c0-01e10e6a3677
# ╠═5042ab46-f698-11ea-3847-4b6867cdaa18
