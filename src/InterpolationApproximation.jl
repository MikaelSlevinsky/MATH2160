const math2160 = joinpath(dirname(@__DIR__), "figures")

import Pkg

Pkg.add("PyPlot")

using PyPlot

function lagrange(knots, data)
#LAGRANGE  Plots the Lagrange polynomial interpolant for the
#given DATA at the given KNOTS

    p = lagrange_poly(knots, data);

    left = minimum(knots); right = maximum(knots); width = right - left;
    top = maximum(data); bottom = minimum(data); height = top - bottom;

    x = range(left - 0.1*width, stop = right + 0.1*width, length = 2000);
    y = polyval(p, x);

    T = Float64
    clf()
    plot(T.(knots), T.(data), "mx")
    plot(T.(x), T.(y))
    xlabel("\$x\$")
    xlim((T(x[1]),T(x[end])))
    ylim((T(bottom-0.1*height),T(top+0.1*height)))
    #axis((T(left-0.1*width),T(right+0.1*width),T(bottom-0.1*height),T(top+0.1*height)))
    axis((T(left),T(right),T(bottom-0.1*height),T(top+0.1*height)))
    p
end

function lagrange_poly(knots, data)
#LAGRANGE_POLY  Lagrange interpolation polynomial.
#   P = LAGRANGE_POLY(KNOTS, DATA) returns a vector containing the
#   coefficients of the Lagrange polynomial for the function
#   values in the vector DATA evaluated at the points KNOTS.
#
#   The N coefficients are ordered so that P(1) is the
#   coefficient corresponding to x^(N-1) and P(N) is the
#   coefficient for x^0, i.e., the polynomial l(x) is defined
#   by l(x) = P(1)x^(N-1) + P(2)x^(N-2) + ... + P(N)x^0.
#
#   The function computes the coefficients of the Lagrange
#   polynomial using the constructive method shown in lectures.

    p = zero(data)
    for k = 1:length(knots)
        p = p + generate_L(k, knots)*data[k]
    end
    p
end

function generate_L(k, knots)
# Subfunction for generating the individual components
    these_points = vcat(knots[1:(k-1)], knots[(k+1):end])
    Lk = poly(these_points)
    the_scale_factor = prod(knots[k] .- these_points)
    Lk = Lk/the_scale_factor
    Lk
end

function poly(z)
    n = length(z)
    p = zeros(eltype(z), n+1)
    p[1] = 1
    for j = 1:n
        p[2:j+1] = p[2:j+1]-z[j]*p[1:j]
    end
    p
end

function polyval(p::Vector,z::Number)
    ret = p[1]
    for k = 2:length(p)
        ret = z*ret + p[k]
    end
    ret
end
polyval(p::Vector,z::AbstractVector) = [ polyval(p,zi) for zi in z ]

function chebyshevpoints(::Type{T}, n::Integer; kind::Integer=1) where T <: Number
    if kind == 1
        T[sinpi((n-2k-one(T))/2n) for k=n-1:-1:0]
    elseif kind == 2
        if n == 1
            zeros(T,1)
        else
            T[sinpi((n-2k-one(T))/(2n-2)) for k=n-1:-1:0]
        end
    end
end
chebyshevpoints(n::Integer; kind::Integer=1) = chebyshevpoints(Float64, n; kind=kind)

f = x-> 1/(1+x^2)

for i=1:3
    local x = -5:big(2.0)/2^i:5
    lagrange(x,f.(x));ylabel("\$f(x)\$");grid(true);gcf()
    savefig(math2160*"/runge$i.pdf")
end

for i=1:3
    local x = 5chebyshevpoints(BigFloat,1+5*2^i)
    lagrange(x,f.(x));ylabel("\$f(x)\$");grid(true);gcf()
    savefig(math2160*"/rungefixed$i.pdf")
end

function chebyshevt(k::Int,x)
    cos(k*acos(x))
end

function chebyshevu(k::Int,x)
    sin((k+1)*acos(x))/sin(acos(x))
end

function legendre(k::Int,x)
    if k == 0
        ret = zero(x) + 1
    elseif k == 1
        ret = x
    else
        t1 = zero(x)+1
        t2 = x
        for i=2:k
            t1, t2 = t2, ((2i-1)*x*t2 - (i-1)*t1)/i
        end
        ret = t2
    end
    ret
end

x = chebyshevpoints(100)
leg = Vector{String}(undef, 6)

clf()
for i=0:5
    plot(x,chebyshevt.(i, x))
    leg[i+1] = "\$T_$(i)(x)\$"
end
legend(leg,loc="lower right");grid(true)
xlabel("\$x\$");ylabel("Chebyshev polynomials of the first kind")
ylim((-1.1,1.1))
savefig(math2160*"/chebyshevt.pdf")

clf()
for i=0:5
    plot(x,chebyshevu.(i, x))
    leg[i+1] = "\$U_$(i)(x)\$"
end
legend(leg,loc="lower right");grid(true)
xlabel("\$x\$");ylabel("Chebyshev polynomials of the second kind")
savefig(math2160*"/chebyshevu.pdf")

clf()
for i=0:5
    plot(x,legendre.(i, x))
    leg[i+1] = "\$P_$(i)(x)\$"
end
legend(leg,loc="lower right");grid(true)
xlabel("\$x\$");ylabel("Legendre polynomials")
ylim((-1.1,1.1))
savefig(math2160*"/legendre.pdf")

#=
using Remez

colours = ["b","g","r","c","m"]
x = chebyshevpoints(501)

f = x -> abs(x+0.5)

clf()
i=1
for n in [1;2;4;8;16]
    (cfsabs,err,equi) = remez(f,(-1,1),n)
    bestapprox = [Float64(Remez.clenshaw(cfsabs, big(x))) for x in x]
    plot(x, bestapprox, colours[i])
    leg[i] = "\$ p_{$(n)}(x)\$"
    i+=1
end
plot(x,f.(x),"-k")
leg[i] = "\$|x+\\frac{1}{2}|\$"
legend(leg,loc="lower right");grid(true)
xlabel("\$x\$");ylabel("Best Polynomial Approximation in \$L^\\infty\\!([-1,1])\$")
savefig(math2160*"/BPAabs.pdf")

clf()
i=1
for n in [1;2;4;8;16]
    (cfsabs,err,equi) = remez(f,(-1,1),n)
    bestapprox = [Float64(Remez.clenshaw(cfsabs, big(x))) for x in x]
    plot(x, f.(x) - bestapprox, colours[i])
    leg[i] = "\$ f(x) - p_{$(n)}(x)\$"
    i+=1
end
legend(leg[1:end-1],loc="lower right");grid(true)
i=1
for n in [1;2;4;8;16]
    (cfsabs,err,equi) = remez(f,(-1,1),n)
    err = Float64(err)
    plot([-1;1],[err;err],"--"*colours[i],[-1;1],[-err;-err],"--"*colours[i])
    i+=1
end
xlabel("\$x\$");ylabel("Best Polynomial Approximation in \$L^\\infty\\!([-1,1])\$")
savefig(math2160*"/BPAerrabs.pdf")
=#

x = -3:3
yi = [ones(4);zeros(3)]

clf()

lagrange(x,yi);xlabel("")
ax = gca()
ax[:set_aspect]("equal")
ax[:spines]["left"][:set_color]("none")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_color]("none")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim(1.02collect(extrema(x)));ylim((-0.3,1.4))
savefig(math2160*"/spline1a.pdf")

clf();PyPlot.axes(aspect="equal")

plot(x,yi,"xm",x,yi)
ax = gca()
ax[:spines]["left"][:set_color]("none")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_color]("none")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim(1.02collect(extrema(x)));ylim((-0.3,1.4))
savefig(math2160*"/spline1b.pdf")


f = x -> x#1/(1+exp(x))

n = 1000
h = 10/n
k = -n:n
kh = k*h
wk = exp.(-kh.^2)
fk = f.(kh)

function cardinalbary(x::Number)
    ind = findfirst(v->x∈v,kh)
    if ind == nothing
        num = (-1.0).^k.*wk.*fk./(x.-kh)
        den = (-1.0).^k.*wk./(x.-kh)
        return sum(num)/sum(den)
    else
        return fk[ind]
    end
end

clf()
x = range(-10, stop = 10, length = 1002)
semilogy(x,abs.(f.(x) .- cardinalbary.(x)),"-k")
grid(true)
xlabel("\$x\$");ylabel("Absolute error")
savefig(math2160*"/cardinalbary.pdf")
