### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 67d2c70a-fe09-11ea-08a3-9304985aad84
using LinearAlgebra

# ╔═╡ 74f1fe0c-fcdc-11ea-008d-13f0aff0f2dd
md"""The Hilbert matrix:

$$H_n = \begin{bmatrix} 1 & \frac{1}{2} & \frac{1}{3} & \cdots & \frac{1}{n}\\ \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \cdots & \frac{1}{n+1}\\ \tfrac{1}{3} & \frac{1}{4} & \frac{1}{5} & \cdots & \frac{1}{n+2}\\ \vdots & \vdots & \vdots & \ddots\\\frac{1}{n} & \frac{1}{n+1} & \frac{1}{n+2} & \cdots & \frac{1}{2n-1}\end{bmatrix}.$$

is a great example to illustrate ill-conditioning. For a symmetric ill-conditioned matrix, it is useful to examine sensitivities round eigenvectors as they all "point" in different directions, colloquially speaking."""

# ╔═╡ e2c4951c-fcdb-11ea-32ed-0166c7d10495
H = n -> inv.((1:n) .+ (1:n)' .- 1)

# ╔═╡ ab6521f4-fcdb-11ea-08a9-7dc726d94453
A = H(15)

# ╔═╡ 747ec76a-fe09-11ea-264d-7f286c10e635
Λ, V = eigen(A)

# ╔═╡ 755ca21a-fe09-11ea-2d4b-c5d66090455e
x = V[:,1]

# ╔═╡ 6e6f031c-fe09-11ea-2f32-456fe3fb5b64
Δx = 1e-6V[:, 8]

# ╔═╡ d2af3752-fe09-11ea-150e-173837f1d9a1
md"""In the relative sense, the perturbation $\Delta x$ is small in the $2$-norm."""

# ╔═╡ ccadaf34-fe09-11ea-0159-33bfd44a6573
norm((x+Δx)-x)/norm(x)

# ╔═╡ f08bb430-fe09-11ea-3b5e-fb9e8577b8c3
md"""But the image of the perturbed $x+\Delta x$ under $H_{15}$ is much farther away from the image of $x$."""

# ╔═╡ ee9b9f0a-fe09-11ea-043e-23819b897ffc
y = A*x

# ╔═╡ 2535cdd8-fe0a-11ea-1337-dd7822f4756a
Δy = A*Δx

# ╔═╡ 2baae568-fe0a-11ea-1039-7348579cb424
norm(y+Δy-y)/norm(y)

# ╔═╡ 32ea7b36-fe0a-11ea-18c1-25e2f8ed78e4
md"""It is always true that this relative perturbation is bounded above by the condition number of $H_{15}$ multiplied by the relative perturbation in the domain."""

# ╔═╡ eb9a9d2e-fe09-11ea-34e6-ab182f4f59aa
norm(y+Δy-y)/norm(y) ≤ cond(A)*norm((x+Δx)-x)/norm(x)

# ╔═╡ 70a04a84-fe0a-11ea-3667-e534fe144c8e
md"""This is because the condition number of the Hilbert matrix is indeed quite large."""

# ╔═╡ 7ea90c86-fe0a-11ea-3d31-9922eca502af
cond(A)

# ╔═╡ Cell order:
# ╟─74f1fe0c-fcdc-11ea-008d-13f0aff0f2dd
# ╠═e2c4951c-fcdb-11ea-32ed-0166c7d10495
# ╠═ab6521f4-fcdb-11ea-08a9-7dc726d94453
# ╠═67d2c70a-fe09-11ea-08a3-9304985aad84
# ╠═747ec76a-fe09-11ea-264d-7f286c10e635
# ╠═755ca21a-fe09-11ea-2d4b-c5d66090455e
# ╠═6e6f031c-fe09-11ea-2f32-456fe3fb5b64
# ╟─d2af3752-fe09-11ea-150e-173837f1d9a1
# ╠═ccadaf34-fe09-11ea-0159-33bfd44a6573
# ╟─f08bb430-fe09-11ea-3b5e-fb9e8577b8c3
# ╠═ee9b9f0a-fe09-11ea-043e-23819b897ffc
# ╠═2535cdd8-fe0a-11ea-1337-dd7822f4756a
# ╠═2baae568-fe0a-11ea-1039-7348579cb424
# ╟─32ea7b36-fe0a-11ea-18c1-25e2f8ed78e4
# ╠═eb9a9d2e-fe09-11ea-34e6-ab182f4f59aa
# ╟─70a04a84-fe0a-11ea-3667-e534fe144c8e
# ╠═7ea90c86-fe0a-11ea-3d31-9922eca502af
