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

# ╔═╡ 520b4d48-f6a7-11ea-3931-53d0163fb2ea
md"""In Julia, there are many different number types. All are `subtypes` of the abstract `supertype` `Number`."""

# ╔═╡ 73dc208c-f6a7-11ea-276a-b7b04280c752
subtypes(Number)

# ╔═╡ 065facb4-f6a8-11ea-0d11-b5c847b9c267
md"""Subtypes can also be abstract, concrete, or parametric."""

# ╔═╡ 2d90f778-f6a8-11ea-1085-1519fa860784
isabstracttype.(subtypes(Number))

# ╔═╡ 63be34fa-f6a8-11ea-0544-cbb6c1096455
isconcretetype.(subtypes(Number))

# ╔═╡ 784b6516-f6a8-11ea-0ec2-6fcd43e7f7f4
md"""We can conclude that `Complex` is a parametric type while `Real` is an abstract type. Usually, but not always, abstract types have subtypes."""

# ╔═╡ 73c0b7fc-f6a7-11ea-258f-8761a263a8c4
subtypes(Real)

# ╔═╡ 738d2ea0-f6a7-11ea-24ac-b56257cafe4c
subtypes(AbstractFloat)

# ╔═╡ 736a6398-f6a7-11ea-1af6-45502f94b8c7
subtypes(Integer)

# ╔═╡ 73433c1e-f6a7-11ea-3e40-1da6858bfa10
subtypes(Signed)

# ╔═╡ be72f4e0-f6a7-11ea-1008-ad6949ca3e2d
subtypes(Unsigned)

# ╔═╡ cdc72fe2-f6a7-11ea-0ef2-c73c54e02c6f
Complex

# ╔═╡ d0b1ed62-f6a7-11ea-06c7-b978a3adf2e2
AbstractIrrational

# ╔═╡ e325ed2e-f6a7-11ea-0acf-f9fafdfb4a0e
Rational

# ╔═╡ 03997c50-f6a9-11ea-0195-1bd12117c044
md"""Mathematically, it's not inconceivable that would wish to work with a rational number type. There are a few problems with this when it comes to arithmetic on a computer. The main issue is that of overflow and underflow. Take, for example, the Hilbert matrix:

$$H_n = \begin{bmatrix} 1 & \frac{1}{2} & \frac{1}{3} & \cdots & \frac{1}{n}\\ \frac{1}{2} & \frac{1}{3} & \frac{1}{4} & \cdots & \frac{1}{n+1}\\ \tfrac{1}{3} & \frac{1}{4} & \frac{1}{5} & \cdots & \frac{1}{n+2}\\ \vdots & \vdots & \vdots & \ddots\\\frac{1}{n} & \frac{1}{n+1} & \frac{1}{n+2} & \cdots & \frac{1}{2n-1}\end{bmatrix}.$$

This matrix is easy enough to create. Julia even allows us to find its inverse with rationals. The catch is that if $n$ is too large, this seemingly innocent-looking matrix's inverse is no longer representable as a ratio of two $64$-bit integers. It could be done with arbitrary precision, but this comes at a significant computational expense.
"""

# ╔═╡ eaf812d4-f6a9-11ea-203e-edfff1f945be
H = n -> inv.((1:n) .+ (1:n)' .- 1)

# ╔═╡ ec1b4bfe-f6a9-11ea-0c9a-63098097eb2f
H(5//1)

# ╔═╡ ecc6e5fc-f6a9-11ea-0e82-1b1287d70741
inv(H(5//1))

# ╔═╡ ef935fea-f6a9-11ea-1e87-63002868c0bb
inv(H(5//1))*H(5//1)

# ╔═╡ ef542708-f6a9-11ea-1373-a789a8d87284
H(15//1)

# ╔═╡ ef3a3c80-f6a9-11ea-2ce4-4d5ba1f1c1b7
inv(H(15//1))

# ╔═╡ d93f45b4-f6aa-11ea-3c46-d3dab68b6bb5
md"""This is one reason we tend to use floating-point types and arithmetic in numerical analysis."""

# ╔═╡ 69204d64-f6aa-11ea-00fb-5d5574b9e0a8
H(15)

# ╔═╡ ef1b3984-f6a9-11ea-1851-05430f4c1766
inv(H(15))

# ╔═╡ efa59b50-f6aa-11ea-1ea5-754c31aaf5d3
bitstring(1.0)

# ╔═╡ cbbc6066-f6aa-11ea-05fd-ef11b11356a2
2^9+2^8+2^7+2^6+2^5+2^4+2^3+2^2+2^1+2^0

# ╔═╡ bf73b6d8-f6aa-11ea-3246-a11e5458d7e0
begin
	function colorbitstring(x::Float16)
		s = bitstring(x)
		HTML("""<potato style="color:red">$(string(s[1]))</potato><potato style="color:green">$(s[2:6])</potato><potato style="color:blue">$(s[7:end])</potato>""")
	end
	function colorbitstring(x::Float32)
		s = bitstring(x)
		HTML("""<potato style="color:red">$(string(s[1]))</potato><potato style="color:green">$(s[2:9])</potato><potato style="color:blue">$(s[10:end])</potato>""")
	end
	function colorbitstring(x::Float64)
		s = bitstring(x)
		HTML("""<potato style="color:red">$(string(s[1]))</potato><potato style="color:green">$(s[2:12])</potato><potato style="color:blue">$(s[13:end])</potato>""")
	end
end

# ╔═╡ e7028138-f6ab-11ea-1fba-5b07169e4b67
x = 1.0

# ╔═╡ bf3d4e7c-f6aa-11ea-0b4b-77934236be8b
colorbitstring(Float16(x))

# ╔═╡ d221ccec-f6ab-11ea-3424-b5cbecdea12d
colorbitstring(Float32(x))

# ╔═╡ d828f6e0-f6ab-11ea-30a3-1514b9fad105
colorbitstring(x)

# ╔═╡ bf081f86-f6aa-11ea-270e-21924da7e375
md"""Floating-point is not perfect; it inherently comes with a so-called rounding or pruning of the least significant information. The most shocking way to see this is with $0.1+0.2$:"""

# ╔═╡ beee6f1c-f6aa-11ea-078b-15357188dd68
0.1+0.2

# ╔═╡ 4ae7a3ae-f6ac-11ea-27d2-0dda73e982b8
@bind b html"""<input type=range min=-53 max=53 value=-53>"""

# ╔═╡ 4acaa90c-f6ac-11ea-1729-bb0592e19764
ϵ = 2.0^b

# ╔═╡ 9b29c3ec-f6ac-11ea-1ea5-2bb886d822b1
1+ϵ

# ╔═╡ 4aaf2c2c-f6ac-11ea-23b9-9fb8872e86c4
colorbitstring(Float16(1.0+ϵ))

# ╔═╡ 4a9d590c-f6ac-11ea-3c96-9d3597e7c055
colorbitstring(Float32(1.0+ϵ))

# ╔═╡ 8adb8426-f6ac-11ea-3263-978408076918
colorbitstring(1.0+ϵ)

# ╔═╡ Cell order:
# ╟─520b4d48-f6a7-11ea-3931-53d0163fb2ea
# ╠═73dc208c-f6a7-11ea-276a-b7b04280c752
# ╟─065facb4-f6a8-11ea-0d11-b5c847b9c267
# ╠═2d90f778-f6a8-11ea-1085-1519fa860784
# ╠═63be34fa-f6a8-11ea-0544-cbb6c1096455
# ╟─784b6516-f6a8-11ea-0ec2-6fcd43e7f7f4
# ╠═73c0b7fc-f6a7-11ea-258f-8761a263a8c4
# ╠═738d2ea0-f6a7-11ea-24ac-b56257cafe4c
# ╠═736a6398-f6a7-11ea-1af6-45502f94b8c7
# ╠═73433c1e-f6a7-11ea-3e40-1da6858bfa10
# ╠═be72f4e0-f6a7-11ea-1008-ad6949ca3e2d
# ╠═cdc72fe2-f6a7-11ea-0ef2-c73c54e02c6f
# ╠═d0b1ed62-f6a7-11ea-06c7-b978a3adf2e2
# ╠═e325ed2e-f6a7-11ea-0acf-f9fafdfb4a0e
# ╟─03997c50-f6a9-11ea-0195-1bd12117c044
# ╠═eaf812d4-f6a9-11ea-203e-edfff1f945be
# ╠═ec1b4bfe-f6a9-11ea-0c9a-63098097eb2f
# ╠═ecc6e5fc-f6a9-11ea-0e82-1b1287d70741
# ╠═ef935fea-f6a9-11ea-1e87-63002868c0bb
# ╠═ef542708-f6a9-11ea-1373-a789a8d87284
# ╠═ef3a3c80-f6a9-11ea-2ce4-4d5ba1f1c1b7
# ╟─d93f45b4-f6aa-11ea-3c46-d3dab68b6bb5
# ╠═69204d64-f6aa-11ea-00fb-5d5574b9e0a8
# ╠═ef1b3984-f6a9-11ea-1851-05430f4c1766
# ╠═efa59b50-f6aa-11ea-1ea5-754c31aaf5d3
# ╠═cbbc6066-f6aa-11ea-05fd-ef11b11356a2
# ╠═bf73b6d8-f6aa-11ea-3246-a11e5458d7e0
# ╠═e7028138-f6ab-11ea-1fba-5b07169e4b67
# ╠═bf3d4e7c-f6aa-11ea-0b4b-77934236be8b
# ╠═d221ccec-f6ab-11ea-3424-b5cbecdea12d
# ╠═d828f6e0-f6ab-11ea-30a3-1514b9fad105
# ╟─bf081f86-f6aa-11ea-270e-21924da7e375
# ╠═beee6f1c-f6aa-11ea-078b-15357188dd68
# ╠═4ae7a3ae-f6ac-11ea-27d2-0dda73e982b8
# ╠═4acaa90c-f6ac-11ea-1729-bb0592e19764
# ╠═9b29c3ec-f6ac-11ea-1ea5-2bb886d822b1
# ╠═4aaf2c2c-f6ac-11ea-23b9-9fb8872e86c4
# ╠═4a9d590c-f6ac-11ea-3c96-9d3597e7c055
# ╠═8adb8426-f6ac-11ea-3263-978408076918
