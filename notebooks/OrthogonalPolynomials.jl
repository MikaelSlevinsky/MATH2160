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

# ╔═╡ 323a15ee-1858-11eb-2b20-6d60ad1c40a9
using ApproxFun, LinearAlgebra, Plots

# ╔═╡ dcdd7774-187c-11eb-02b6-096b31338fa1
md"""[ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl) is a Julia package for numerically computing with functions. The core concept is that it represents a function by a Julia struct `Fun` with two fields: `coefficients` to store the floating-point data; and `space` to distinguish how the coefficients are interpreted.

Let's start by creating the identity function."""

# ╔═╡ a9684872-185b-11eb-0d21-3b1ba212ef85
x = Fun(identity)

# ╔═╡ bc09e2c2-185b-11eb-004a-f30cf0d3f94a
plot(x; label="Id(x)")

# ╔═╡ 3f926d0c-187d-11eb-0d5b-a75a1a3c98e7
md"""The default behaviour is to use the `Chebyshev()` approximation space so that the coefficients are first kind Chebyshev coefficients of the function `x`. A `Fun` `f` in the `Chebyshev()` space is the numerical implementation of:

$$f(x) = \sum_{k=0}^n c_k T_k(x)\qquad T_k(x) = \cos(k\cos^{-1}(x)).$$

The `Fun` can be evaluated at any point in its domain."""

# ╔═╡ 3f7306a6-187d-11eb-29b9-317c27108cb3
domain(x)

# ╔═╡ 3f57fbae-187d-11eb-3996-57d15caac609
x(0.5)

# ╔═╡ 5e73ee50-1880-11eb-37d0-336ea444dde5
md"""or extrapolation can be used outside $[-1,1]$."""

# ╔═╡ 6b731f72-1880-11eb-0df2-0f30dfab433f
extrapolate(x, 1.5)

# ╔═╡ 3f353c0e-187d-11eb-3f0f-d3bc12033168
md"""It can also be differentiated by using the apostrophe, `'`."""

# ╔═╡ 3f1a8c8a-187d-11eb-0853-cf4a184e528d
x'

# ╔═╡ 3ef7d314-187d-11eb-38e0-7936b0e44ca2
md"""Let's now use `x` to create a weight function for which we will compute the first few orthonormal polynomials through the function `lanczos`, which implements Gram--Schmidt orthonormalization (compare Lemma 3.2.6 to the [source](https://github.com/JuliaApproximation/ApproxFun.jl/blob/master/src/Extras/lanczos.jl))."""

# ╔═╡ f0bad210-185b-11eb-1e79-5baced80673a
w = one(x)

# ╔═╡ fc169e0c-185b-11eb-18f3-a14064c6a773
P, β, γ = lanczos(w, 10)

# ╔═╡ fdb88aaa-187d-11eb-26f1-25fdba7bc61d
md"""We can check that they're degree-graded by using `ncoefficients` to return the number of coefficients in each expansion."""

# ╔═╡ fa90bd84-187d-11eb-3174-01684dc378bf
[ncoefficients(P[j]) for j in 1:length(P)]

# ╔═╡ b3456128-187d-11eb-224a-33f657eb06ab
md"""We can also check that they're orthonormal by using `sum`, which has a method to mean "definite integral." We'll compare the array of inner products against the identity and take the Frobenius norm of the difference."""

# ╔═╡ 1cf8c336-185c-11eb-0e75-d9fa1adb37d4
norm([sum(P[k]*P[j]*w) for k in 1:length(P), j in 1:length(P)] - I)

# ╔═╡ f8b1dab6-187d-11eb-10a8-b7a3b93ba1a8
md"""Last but not least, we can plot them, where it is easy to observe the consequence of Theorem 3.2.8: the roots are all in $(-1,1)$. With $w(x) = 1$, these polynomials are the orthonormalized Legendre polynomials, though they are expanded in the first kind Chebyshev polynomial basis."""

# ╔═╡ 5978a2e8-185c-11eb-12cd-875ed64b881b
begin
	p = plot(P[1]; legend=false)
	for k in 2:length(P)
		plot!(pad(P[k], 3*ncoefficients(P[k])))
	end
	p
end

# ╔═╡ b04de64e-187e-11eb-1f04-c1d43a5370fc
md"""Binding a Julia variable `a` to an html widget gives us some freedom to explore orthogonal polynomials reactively. In `ApproxFun`, the weight function can't be just anything, so let me show you a couple examples to give you an idea."""

# ╔═╡ 2b4b458c-185d-11eb-32ec-0b0180f718ce
@bind a html"""<input type='range' min = 0.0 max = 0.9 step = 0.01>"""

# ╔═╡ cd3fd346-187e-11eb-0489-17fc7ef09394
w1 = (1-x)^a

# ╔═╡ cd6e516a-187e-11eb-1798-5b57bdd16fb6
P1, β1, γ1 = lanczos(w1, 10);

# ╔═╡ cd8d0ecc-187e-11eb-3719-670578126046
begin
	p1 = plot(P1[1]; label="\$P_n^{($a, 0)}(x)\$") # Jacobi polynomials
	for k in 2:length(P1)
		plot!(pad(P1[k], 3*ncoefficients(P1[k])); label="")
	end
	p1
end

# ╔═╡ cda6caba-187e-11eb-313c-b14ddc46c709
w2 = exp(5*a*x) # A non-classical weight.

# ╔═╡ cdc09ce2-187e-11eb-3842-c90047c897f0
P2, β2, γ2 = lanczos(w2, 10);

# ╔═╡ cdda45c2-187e-11eb-327c-ab884423e7b7
begin
	p2 = plot(P2[1]; legend=false)
	for k in 2:length(P2)
		plot!(pad(P2[k], 3*ncoefficients(P2[k])))
	end
	p2
end

# ╔═╡ cdf092a8-187e-11eb-0cb2-9517fe650af4
w3 = Fun(one, Chebyshev(-1..(-a))∪Chebyshev(a..1)) # A constant weight on a piecewise-defined domain, also non-classical.

# ╔═╡ e3db2938-187f-11eb-32e6-9d97ebe3a117
P3, β3, γ3 = lanczos(w3, 10);

# ╔═╡ ea48f8cc-187f-11eb-1d5f-fb02f9d576a1
begin
	p3 = plot(P3[1]; legend=false)
	for k in 2:length(P3)
		plot!(pad(P3[k], 3*ncoefficients(P3[k])))
	end
	p3
end

# ╔═╡ Cell order:
# ╠═323a15ee-1858-11eb-2b20-6d60ad1c40a9
# ╟─dcdd7774-187c-11eb-02b6-096b31338fa1
# ╠═a9684872-185b-11eb-0d21-3b1ba212ef85
# ╠═bc09e2c2-185b-11eb-004a-f30cf0d3f94a
# ╟─3f926d0c-187d-11eb-0d5b-a75a1a3c98e7
# ╠═3f7306a6-187d-11eb-29b9-317c27108cb3
# ╠═3f57fbae-187d-11eb-3996-57d15caac609
# ╟─5e73ee50-1880-11eb-37d0-336ea444dde5
# ╠═6b731f72-1880-11eb-0df2-0f30dfab433f
# ╟─3f353c0e-187d-11eb-3f0f-d3bc12033168
# ╠═3f1a8c8a-187d-11eb-0853-cf4a184e528d
# ╟─3ef7d314-187d-11eb-38e0-7936b0e44ca2
# ╠═f0bad210-185b-11eb-1e79-5baced80673a
# ╠═fc169e0c-185b-11eb-18f3-a14064c6a773
# ╟─fdb88aaa-187d-11eb-26f1-25fdba7bc61d
# ╠═fa90bd84-187d-11eb-3174-01684dc378bf
# ╟─b3456128-187d-11eb-224a-33f657eb06ab
# ╠═1cf8c336-185c-11eb-0e75-d9fa1adb37d4
# ╟─f8b1dab6-187d-11eb-10a8-b7a3b93ba1a8
# ╠═5978a2e8-185c-11eb-12cd-875ed64b881b
# ╟─b04de64e-187e-11eb-1f04-c1d43a5370fc
# ╠═2b4b458c-185d-11eb-32ec-0b0180f718ce
# ╠═cd3fd346-187e-11eb-0489-17fc7ef09394
# ╠═cd6e516a-187e-11eb-1798-5b57bdd16fb6
# ╠═cd8d0ecc-187e-11eb-3719-670578126046
# ╠═cda6caba-187e-11eb-313c-b14ddc46c709
# ╠═cdc09ce2-187e-11eb-3842-c90047c897f0
# ╠═cdda45c2-187e-11eb-327c-ab884423e7b7
# ╠═cdf092a8-187e-11eb-0cb2-9517fe650af4
# ╠═e3db2938-187f-11eb-32e6-9d97ebe3a117
# ╠═ea48f8cc-187f-11eb-1d5f-fb02f9d576a1
