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

# ╔═╡ af503bac-f306-11ea-093a-75c6cf4beaf5
using Plots

# ╔═╡ cc327ac8-f306-11ea-03a9-01767dcdf51b
using LinearAlgebra

# ╔═╡ bd7c0d78-f306-11ea-0967-e329242407af
data = [1,0,1,1,0,1,2,1,5,4,3,3,4,3,17,6,6,11,18,15,32,56,54,89,100,157,129,146,214,241,142,621,701,617,634,714,898,665,1128,1164,1119,1552,1092,1537,1600,1155,1230,1541,1327,1383,1170,1065,1297,1383,1316,1727,1821,1456,1673,1773,1593,1768,1920,1778,1466,1541,1605,1526,1571,1639,1825,1653,2760,1298,1274,1450,1426,1512,1268,1146,1133,1176,1121,1123,1212,1251,1138,1070,1040,1030,1182,1156,1141,1078,1012,936,872,993,906,772,757,758,705,675,641,609,722,642,545,409,472,405,413,467,377,360,320,386,367,409,390,318,300,326,279,380,172,238,218,668,286,67,501,319,226,219,399,232,267,371,321,221,243,565,331,343,435,405,330,339,786,573,543,432,534,350,355,686,397,476,329,513,287,285,147,761,395,374,424,236,230,681,289,423,390,418,237,198,785,282,336,383,499,257,267,751,322,448,431,510,315,267,1008,477,498,570,631,371,400,247,1606,546]

# ╔═╡ 25ad7a48-f308-11ea-3e74-b9f58dfa89cd
N = length(data)

# ╔═╡ cccd1574-f306-11ea-1503-a12916d3c16c
begin
	days = 1:N
	scaled_days = -1 .+ 2/(N-1).*(days.-1)
end

# ╔═╡ 2f7fbba4-f307-11ea-179f-e5a750fbfcd9
@bind deg html"""<input type="range" min="0" max="201" value=5>"""

# ╔═╡ 462ba9ee-f307-11ea-00d7-0ba8d55c2ea9
deg

# ╔═╡ cc9b3f36-f306-11ea-38ed-496a1dd0356b
begin
	A = scaled_days.^(0:deg)'
	c = A\data
end

# ╔═╡ cc81c696-f306-11ea-371f-3bce5f661229
f = x-> c[1] + x*(c[2] + x*(c[3] + x*(c[4] + x*(c[5] + x*c[6])))) # A degree-5 polynomial with coefficients `c` evaluated at `x`.

# ╔═╡ cc65fb96-f306-11ea-0075-d14d27f13b3b
function horner(x, c)
	N = length(c)
	ret = c[N]
	for k in N-1:-1:1
		ret = x*ret+c[k]
	end
	return ret
end

# ╔═╡ ccb6c38c-f306-11ea-3782-bf4bf252efd8
begin
	scatter(days, data;label="Real Dataᵀᴹ")
	ylims!(extrema(data))
	plot!(days, map(x->horner(x,c), scaled_days); label="Polynomial approximation")
	xlabel!("Days since first case")
	ylabel!("Number of cases")
	title!("Daily New Cases in Canada")
end

# ╔═╡ cc4c4aa2-f306-11ea-04dc-0314fefcf7fe
md" Linear algebra is all about $Ax = b$.

Calculus is all about $f(x)$, $f'(x)$ and $\int_a^bf(x){\rm\,d}x$"

# ╔═╡ a3878152-f30a-11ea-119c-7508a00b83f1
md"For least-squares problems, the norm of the residual, $r = b-Ax$, is not usually small. But if we multiply it by $A^*$, then it is near machine precision multiplied by the condition number of the matrix and the norm of the right-hand side."		

# ╔═╡ ab1d9616-f30a-11ea-0e22-1bd11e77e625
norm(A*c-data) # Pretty huge

# ╔═╡ abb15fa4-f30a-11ea-086c-7127ae89c74c
norm(A'*(A*c-data)) ≤ 2*eps()*cond(A)*norm(data)

# ╔═╡ Cell order:
# ╠═af503bac-f306-11ea-093a-75c6cf4beaf5
# ╠═bd7c0d78-f306-11ea-0967-e329242407af
# ╠═25ad7a48-f308-11ea-3e74-b9f58dfa89cd
# ╠═cccd1574-f306-11ea-1503-a12916d3c16c
# ╠═ccb6c38c-f306-11ea-3782-bf4bf252efd8
# ╠═2f7fbba4-f307-11ea-179f-e5a750fbfcd9
# ╠═462ba9ee-f307-11ea-00d7-0ba8d55c2ea9
# ╠═cc9b3f36-f306-11ea-38ed-496a1dd0356b
# ╠═cc81c696-f306-11ea-371f-3bce5f661229
# ╠═cc65fb96-f306-11ea-0075-d14d27f13b3b
# ╟─cc4c4aa2-f306-11ea-04dc-0314fefcf7fe
# ╠═cc327ac8-f306-11ea-03a9-01767dcdf51b
# ╟─a3878152-f30a-11ea-119c-7508a00b83f1
# ╠═ab1d9616-f30a-11ea-0e22-1bd11e77e625
# ╠═abb15fa4-f30a-11ea-086c-7127ae89c74c
