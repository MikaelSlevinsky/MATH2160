const math2160 = "/Users/Mikael/Documents/Courses/MATH2160"

using PyPlot

f = x -> exp(x)
fp = x -> exp(x)

x = 1.0
h = logspace(log10(eps()/2),0,1000)#2.^(0:-1.0:-52)

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

a,b,h = -1.0,1.0,logspace(-7,0,250)
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


using DualNumbers

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

function divideddiff(f::Vector, x::Vector)
    @assert length(f) == length(x)
    T, n = promote_type(eltype(f), eltype(x)), length(x)
    df = zeros(T, n)
    for i=1:n
        df[i] = f[i]
    end
    for j = 2:n, i = n:-1:j
        df[i] = (df[i] - df[i-1])/(x[i] - x[i-j+1])
    end
    df
end

richardson(a::Vector, h::Vector) = divideddiff(a ./ h, h)./divideddiff(1 ./ h, h)

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
a = 1./n.^2
c = cumsum(a)
h = 1./n
richardson(c[1:5:51],h[1:5:51])-π^2/6


function forwardeuler(f,y0,n,a,b)
    h = (b-a)/(n-1)
    t = collect(linspace(a,b,n))
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
    t = collect(linspace(a,b,n))
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
h = 2.0.^(-n)
for i in n
    retf[i] = forwardeuler((t,y) -> y,1.0,2^i+1,0,1)[2][end]
    retb[i] = backwardeuler(1.0,2^i+1,0,1)[2][end]
end

richardson(retf,h) - e
richardson(retb,h) - e


# Gauss--Legendre
function gausslegendre(n)
    sqrtbeta = sqrt.((1:n-1).^2 ./ (4 .* (1:n-1) .^ 2 .- 1))
    μ0 = 2.0
    J = diagm(-1 => sqrtbeta, 1 => sqrtbeta)
    x, Q = eigen(J)
    w = vec(μ0*Q[1,:].^2)
    x, w
end

# Gauss--Hermite
function gausshermite(n)
    sqrtbeta = sqrt.((1:n-1)/2)
    μ0 = sqrt(pi)
    J = diagm(-1 => sqrtbeta, 1 => sqrtbeta)
    x, Q = eigen(J)
    w = vec(μ0*Q[1,:].^2)
    x, w
end

gausshermite(10)


function montecarlo(f,d,N)
    x = [rand(d) for i=1:N]
    val = f(x[1])
    for i=2:N
        val += f(x[i])
    end
    val/N
end

immutable BoxIntegrand{T,d}
    s::T
    BoxIntegrand{T,d}(s::T) where {T,d} = new{T,d}(s)
end

function (B::BoxIntegrand{T,d}){T,d}(x::Matrix{T},i::Int)
    ret = zero(T)
    for j=1:size(x,1)
        ret += x[j,i]^2
    end
    ret^(B.s/2)
end

function montecarlo{d}(f::BoxIntegrand{Float64,d},n::Int)
    x = rand(d,n)
    ret = f(x,1)
    for i=2:n
        ret += f(x,i)
    end
    return ret/n
end

function montecarloantithetic{d}(f::BoxIntegrand{Float64,d},n::Int)
    x = rand(d,n)
    ret = f(x,1)
    for i=2:n
        ret += f(x,i)
    end
    y = 1-x
    for i=1:n
        ret += f(y,i)
    end
    return (ret/2n)
end

B4n2 = π*log(2+sqrt(3))-2catalan-π^2/8
asy = n -> sqrt(n/3)*(1-1/10n)


Bi = BoxIntegrand{Float64,4}(-2.0)
MC = Float64[montecarlo(Bi,i)-B4n2 for i in 10.^(1:7)]
MCA = Float64[montecarloantithetic(Bi,i)-B4n2 for i in 10.^(1:7)]

@show [MC MCA]

for n in [10,1000,1000_000]
    Bi = BoxIntegrand{Float64,n}(1.0)
    MC = Float64[montecarlo(Bi,i)-asy(n) for i in 10.^(1:(8-round.(Int,log10(n))))]
    @show MC
end


function generalizedhorner(xpts,fxpts,x)
    N = length(fxpts)

    ret = fxpts[N]*one(x)
    for i = N-1:-1:1
        ret = (x-xpts[i])*ret + fxpts[i]
    end

    ret
end
