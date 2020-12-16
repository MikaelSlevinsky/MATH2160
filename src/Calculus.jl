const math2160 = joinpath(dirname(@__DIR__), "figures")

import Pkg

Pkg.add(["DualNumbers", "FastGaussQuadrature", "PyPlot"])

using DualNumbers, FastGaussQuadrature, LinearAlgebra, PyPlot

f = x -> exp(x)
fp = x -> exp(x)

x = 1.0
h = exp10.(range(log10(eps()/2), stop = 0, length=1000))

FD1 = [abs(fp(x) - (f(x+h[i])-f(x))/h[i]) for i in eachindex(h)]
FD2 = [abs(fp(x) - (f(x+h[i])-f(x-h[i]))/2h[i]) for i in eachindex(h)]
CD2 = [abs(fp(x) - imag(f(x+im*h[i]))/h[i]) for i in eachindex(h)]

firstorder = h
secondorder = h.^2

clf()
loglog(h,FD1,"-r",h,FD2,"-g",h,CD2,"-b",h,firstorder,"--k",h,secondorder,"--k")
legend(["First-order","Second-order","Complex second-order"],loc="lower left")
xlabel("\$h\$");ylabel("Absolute Error")
xlim(extrema(h));ylim((1e-16,1e1))
annotate("\$\\mathcal{O}(h)\$",[1e-6,1e-4],fontsize=18)
annotate("\$\\mathcal{O}(h^2)\$",[1e-4,1e-10],fontsize=18)
grid(true);gcf()
savefig(math2160*"/finitedifferences.pdf")

function trap(f,a,b,n)
    s = (f(a)+f(b))/2
    h = (b-a)/n
    for i=1:n-1
        s += f(a+i*h)
    end
    h*s
end

function simpson(f,a,b,n)
    @assert iseven(n)
    s = f(a)+f(b)
    h = (b-a)/n
    for i=1:2:n
        s += 4f(a+i*h)
    end
    for i=2:2:n-1
        s += 2f(a+i*h)
    end
    (h/3)*s
end

#n = 2.^(1:25)

a,b,h = -1.0,1.0,exp10.(range(-7, stop = 0, length = 250))
n = 2round.([Int], (b-a)./2h)

CT = [abs(2sinh(1) - trap(f,a,b,i)) for i in n]
CS = [abs(2sinh(1) - simpson(f,a,b,i)) for i in n]

h = (b-a)./n
secondorder = h.^2/2
fourthorder = h.^4/20

clf()
loglog(h,CT,"-r",h,CS,"-g",h,secondorder,"--k",h,fourthorder,"--k")
legend(["Composite Trapezoidal","Composite Simpson",],loc="upper left")
xlabel("\$h\$");ylabel("Absolute Error")
xlim(extrema(h));ylim((1e-16,1))
annotate("\$\\mathcal{O}(h^2)\$",[5e-4,1e-4],fontsize=18)
annotate("\$\\mathcal{O}(h^4)\$",[5e-4,1e-10],fontsize=18)
grid(true);gcf()
savefig(math2160*"/compositenewtoncotes.pdf")


f = x -> sqrt(sin(exp(x))+cos(x^2))
println("This is f(1+ɛ): ",f(1+ɛ))

include("duck.jl")

rsq(t) = xd(t).^2+yd(t).^2
#268779.9125305255

[trap(rsq,0,2π,2^i)/2 for i in 1:10]
[simpson(rsq,0,2π,2^i)/2 for i in 1:10]

rsqBF(t) = xdBF(t).^2+ydBF(t).^2
#268779.9125305254578448085959263059515774569443289752564162584708603294269296908

[trap(rsqBF,big(0),2big(π),2^i)/2 for i in 1:10]
[simpson(rsqBF,big(0),2big(π),2^i)/2 for i in 1:10]

function divideddiff(f::AbstractVector, x::AbstractVector)
    @assert length(f) == length(x)
    T, n = promote_type(eltype(f), eltype(x)), length(x)
    df = zeros(T, n)
    for i = 1:n
        df[i] = f[i]
    end
    for j = 2:n, i = n:-1:j
        df[i] = (df[i] - df[i-1])/(x[i] - x[i-j+1])
    end
    df
end

richardson(a::AbstractVector, h::AbstractVector) = divideddiff(a ./ h, h) ./ divideddiff(1 ./ h, h)

function computecn(::Type{T}, n::Int) where T
    c = zeros(T, n)
    c[1] = 2
    sqrttwo = sqrt(T(2))
    twoi = T(1)
    for i=1:n-1
        twoi *= 2
        c[i+1] = c[i]sqrttwo/sqrt(1+sqrt(1-(c[i]/twoi)^2))
    end
    c
end

c = computecn(Float64, 16)
h = 1 ./ (2 .^ (2:length(c)+1)) .^ 2
rex = richardson(c,h)
for n in eachindex(c)
    println("c[$(lpad(2^n,5))]  =  $(rpad(c[n],18))  and the extrapolant  =  $(rpad(rex[n],18))")
end


n = 1:51
a = inv.(n.^2)
c = cumsum(a)
h = inv.(n)
richardson(c[1:5:51], h[1:5:51]) .- π^2/6

n = 1:10
h1 = inv.(n)
a1 = (1 .+ h1) .^ n
r1 = richardson(a1, h1)
h2 = inv.(2 .^ n)
a2 = (1 .+ h2) .^ (2 .^ n)
r2 = richardson(a2, h2)
l = 19
for n in n
    str = "$(rpad(h1[n], l)) & $(rpad(a1[n], l)) & $(rpad(r1[n], l)) & "
    str *= "$(rpad(h2[n], l)) & $(rpad(a2[n], l)) & $(rpad(r2[n], l))\\\\"
    println(str)
end

n = 1:10
h3 = n.^2 .- 5
a3 = n
r3 = richardson(a3, h3)
l = 19
for n in n
    println("$(rpad(h3[n], 3)) & $(rpad(a3[n], 3)) & $(rpad(r3[n], l))\\\\")
end

y = [2,3,5,7,11,13,17,19,23,29]
z = [5,7,12,16,25,30,39,43,52,65]
x = z.^2 .- 5 .* y.^2
h4 = x./y.^2
a4 = z./y
r4 = richardson(a4, h4)
l = 19
for n in n
    str = "$(rpad(y[n], 3)) & $(rpad(5*y[n]^2, 5)) & $(rpad(z[n]^2, 5)) & "
    str *= "$(rpad(z[n], 3)) & $(rpad(x[n], 3)) & $(rpad(h4[n], l)) & "
    str *= "$(rpad(a4[n], l)) & $(rpad(r4[n], l))\\\\"
    println(str)
end


function forwardeuler(f,y0,n,a,b)
    h = (b-a)/(n-1)
    t = collect(range(a, stop = b, length = n))
    y = zeros(n)
    y[1] = y0
    for i = 1:n-1
        y[i+1] = y[i] + h*f(t[i],y[i])
    end
    t, y
end

# Only for y' = y
function backwardeuler(y0,n,a,b)
    h = (b-a)/(n-1)
    t = collect(range(a, stop = b, length = n))
    y = zeros(n)
    y[1] = y0
    for i = 1:n-1
        y[i+1] = y[i]/(1-h)
    end
    t, y
end

n = 1:15
retf = zeros(length(n))
retb = zeros(length(n))
h = 2.0 .^ -n
for i in n
    retf[i] = forwardeuler((t,y) -> y,1.0,2^i+1,0,1)[2][end]
    retb[i] = backwardeuler(1.0,2^i+1,0,1)[2][end]
end

richardson(retf,h) .- MathConstants.e
richardson(retb,h) .- MathConstants.e


# Gauss--Legendre
function mygausslegendre(n)
    sqrtbeta = sqrt.((1:n-1).^2 ./ (4 .* (1:n-1) .^ 2 .- 1))
    μ0 = 2.0
    J = SymTridiagonal(zeros(n), sqrtbeta)
    x, Q = eigen(J)
    w = vec(μ0*Q[1,:].^2)
    x, w
end

x1, w1 = mygausslegendre(10)
x, w = gausslegendre(10) # from FastGaussQuadrature
hypot(norm(x-x1), norm(w-w1)) < hypot(norm(x), norm(w))

# Gauss--Hermite
function mygausshermite(n)
    sqrtbeta = sqrt.((1:n-1)/2)
    μ0 = sqrt(pi)
    J = SymTridiagonal(zeros(n), sqrtbeta)
    x, Q = eigen(J)
    w = vec(μ0*Q[1,:].^2)
    x, w
end

x1, w1 = mygausshermite(10)
x, w = gausshermite(10) # from FastGaussQuadrature
hypot(norm(x-x1), norm(w-w1)) < hypot(norm(x), norm(w))


l = 19
for n in 1:12
    x, w = gausslegendre(n)
    str = "n = $(lpad(n, 2)), and the integral is approximately: "
    str *= "$(rpad(sum(w.*exp.(-x.^2)), l))"
    println(str)
end
for n in 1:12
    x, w = gausschebyshev(n)
    str = "n = $(lpad(n, 2)), and the integral is approximately: "
    str *= "$(rpad(sum(w./sqrt.(5 .- x.^2)), l))"
    println(str)
end
for n in 1:12
    x, w = gaussjacobi(n, 1/3, 1/3)
    str = "n = $(lpad(n, 2)), and the integral is approximately: "
    str *= "$(rpad(sum(w.*cbrt.(3 .- x))/4, l))"
    println(str)
end
for n in 2 .^ (4:24)
    x, w = gausslegendre(n)
    str = "n = $(lpad(n, 8)), and the integral is approximately: "
    str *= "$(rpad(sum(w.*(sin.(2020*11/2 .* (x .+ 1)) .^ 27))*11/2, l))"
    println(str)
end
for n in 1:20:201
    x, w = gausslaguerre(n)
    str = "n = $(lpad(n, 3)), and the integral is approximately: "
    str *= "$(rpad(sum(w./sqrt.(1 .+ x)), l))"
    println(str)
end
for n in 1:10:81
    x, w = gausslaguerre(n, 1/2)
    f = sqrt.(expm1.(log1p.(x) .* 2 ./ 3)./x)./cbrt.(1 .+ x)
    str = "n = $(lpad(n, 2)), and the integral is approximately: "
    str *= "$(rpad(2/3/exp(1)*sum(w.*f), l))"
    println(str)
end
function hermite(n, x)
    if n == 0
        return one(x)
    elseif n == 1
        return 2*x
    else
        Hkm1 = one(x)
        Hk = 2*x
        for k in 1:n-1
            Hk, Hkm1 = 2*x*Hk-2*k*Hkm1, Hk
        end
        return Hk
    end
end
for n in 9:11
    x, w = gausshermite(n)
    H5 = hermite.(5, x)
    H6 = hermite.(6, x)
    H7 = hermite.(7, x)
    str = "n = $(lpad(n, 2)), and the integral is approximately: "
    str *= "$(rpad(sum(w.*H5.*H6.*H7), l))"
    println(str)
end
for n in 1:25:151
    x, w = gausshermite(n)
    str = "n = $(lpad(n, 3)), and the integral is approximately: "
    str *= "$(rpad(sum(w.*exp.(2 .* x)./(1 .+ x.^2)), l))"
    println(str)
end


function montecarlo(f,d,N)
    x = [rand(d) for i=1:N]
    val = f(x[1])
    for i=2:N
        val += f(x[i])
    end
    val/N
end

struct BoxIntegrand{T, d}
    s::T
    BoxIntegrand{T, d}(s::T) where {T, d} = new{T, d}(s)
end

function (B::BoxIntegrand{T, d})(x::Matrix{T}, i::Int) where {T, d}
    ret = zero(T)
    for j=1:size(x,1)
        ret += x[j,i]^2
    end
    ret^(B.s/2)
end

function montecarlo(f::BoxIntegrand{Float64, d}, n::Int) where d
    x = rand(d,n)
    ret = f(x,1)
    for i=2:n
        ret += f(x, i)
    end
    return ret/n
end

function montecarloantithetic(f::BoxIntegrand{Float64,d}, n::Int) where d
    x = rand(d,n)
    ret = f(x,1)
    for i=2:n
        ret += f(x, i)
    end
    y = 1 .- x
    for i=1:n
        ret += f(y, i)
    end
    return (ret/2n)
end

B4n2 = π*log(2+sqrt(3))-2MathConstants.catalan-π^2/8
asy = n -> sqrt(n/3)*(1-1/10n)

let n = 4
    Bi = BoxIntegrand{Float64,4}(-2.0)
    MC = Float64[montecarlo(Bi,i)-B4n2 for i in 10 .^ (1:7)]
    MCA = Float64[montecarloantithetic(Bi,i)-B4n2 for i in 10 .^ (1:7)]
    @show [MC MCA]
end

for n in [10,1000,1000_000]
    Bi = BoxIntegrand{Float64,n}(1.0)
    MC = Float64[montecarlo(Bi,i)-asy(n) for i in 10 .^ (1:(8-round.(Int, log10(n))))]
    @show MC
end


function generalizedhorner(xpts, fxpts, x)
    N = length(fxpts)
    ret = fxpts[N]*one(x)
    for i = N-1:-1:1
        ret = (x-xpts[i])*ret + fxpts[i]
    end
    ret
end
