const math2160 = "/Users/Mikael/Documents/Courses/MATH2160"

using PyPlot

gray = (0.75,0.75,0.75)
yl = 1.1; xl = golden*yl

ls = cospi.(linspace(0,1,101))

# We first draw the vector norms

clf();axes(aspect="equal")

for (x,y) in (([-1.0;0.0],[0.0;1.0]),([0.0;1.0],[1.0;0.0]),([-1.0;0.0],[0.0;-1.0]),([0.0;1.0],[-1.0;0.0]))
    fill_between(x,y,color=gray)
    plot(x,y,"-k",linewidth=2.0)
end

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norm1.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in ((collect(ls),sqrt.(1-ls.^2)),(collect(ls),-sqrt.(1-ls.^2)))
    fill_between(x,y,color=gray)
    plot(x,y,"-k",linewidth=2.0)
end

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norm2.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in ((collect(ls),(1-ls.^4).^(1/4)),(collect(ls),-(1-ls.^4).^(1/4)))
    fill_between(x,y,color=gray)
    plot(x,y,"-k",linewidth=2.0)
end

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norm4.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in (([-1.0;1.0],[1.0;1.0]),([-1.0;1.0],[-1.0;-1.0]))
    fill_between(x,y,color=gray)
    plot(x,y,"-k",linewidth=2.0)
end
plot([-1.0;-1.0],[-1.0;1.0],"-k",linewidth=2.0)
plot([1.0;1.0],[-1.0;1.0],"-k",linewidth=2.0)

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norminf.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in ((collect(1.6ls),sqrt.(1-ls.^2)),(collect(1.6ls),-sqrt.(1-ls.^2)))
    fill_between(x,y,color=gray)
    plot(x,y,"-k",linewidth=2.0)
end

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/normW.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in (([-1.0;0.0],[0.0;1.0]),([0.0;1.0],[1.0;0.0]),([-1.0;0.0],[0.0;-1.0]),([0.0;1.0],[-1.0;0.0]))
    plot(x,y,"-k",linewidth=2.0)
end
plot([0.0],[1.0],"sk",linewidth=6.0)
plot([1.0],[0.0],"ok",linewidth=6.0)
annotate("\$(0,1)^\\top\$",[0.05,1.05],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.05,0.05],fontsize=18)


ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norm1annotated.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in (([-1.0;2.0],[0.0;2.0]),([1.0;2.0],[0.0;2.0]),([-2.0;-1.0],[-2.0;0.0]),([-2.0;1.0],[-2.0;0.0]))
    plot(x,y,"-k",linewidth=2.0)
end
plot([2.0],[2.0],"sk",linewidth=6.0)
plot([1.0],[0.0],"ok",linewidth=6.0)
annotate("\$(2,2)^\\top\$",[2.1,2.0],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.15,0.1],fontsize=18)


ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-2xl,2xl));ylim((-2yl,2yl))
savefig(math2160*"/normA1annotated.pdf";bbox_inches="tight")


A = [1 2;0 2]
U,Σ,V = svd(A)


clf();axes(aspect="equal")

for (x,y) in ((collect(ls),sqrt.(1-ls.^2)),(collect(ls),-sqrt.(1-ls.^2)))
    plot(x,y,"-k",linewidth=2.0)
end
plot([0.0],[1.0],"sk",linewidth=6.0)
plot([1.0],[0.0],"ok",linewidth=6.0)
plot([0;V[1,1]],[0;V[2,1]],"--k",linewidth=2.0)

annotate("\$(0,1)^\\top\$",[0.05,1.05],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.05,0.05],fontsize=18)

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norm2annotated.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

θ = linspace(0,2,100)
x = Σ[1]*cospi.(θ)
y = Σ[2]*sinpi.(θ)

φ = acos(dot(U[:,1],[1.0;0.0])/norm(U[:,1])/norm([1.0;0.0]))

plot(x*cos(φ)-y*sin(φ),x*sin(φ)+y*cos(φ),"-k",linewidth=2.0)
plot([0;Σ[1]*U[1,1]],[0;Σ[1]*U[2,1]],"--k",linewidth=2.0)

plot([2.0],[2.0],"sk",linewidth=6.0)
plot([1.0],[0.0],"ok",linewidth=6.0)
annotate("\$(2,2)^\\top\$",[2.1,2.05],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.25,0.05],fontsize=18)


ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-2xl,2xl));ylim((-2yl,2yl))
savefig(math2160*"/normA2annotated.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in (([-1.0;1.0],[1.0;1.0]),([-1.0;1.0],[-1.0;-1.0]))
    plot(x,y,"-k",linewidth=2.0)
end
plot([-1.0;-1.0],[-1.0;1.0],"-k",linewidth=2.0)
plot([1.0;1.0],[-1.0;1.0],"-k",linewidth=2.0)
annotate("\$(0,1)^\\top\$",[0.05,1.05],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.05,0.05],fontsize=18)

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-xl,xl));ylim((-yl,yl))
savefig(math2160*"/norminfannotated.pdf";bbox_inches="tight")


clf();axes(aspect="equal")

for (x,y) in (([1.0;3.0],[2.0;2.0]),([-1.0;3.0],[-2.0;2.0]),([-3.0;-1.0],[-2.0;-2.0]),([-3.0;1.0],[-2.0;2.0]))
    plot(x,y,"-k",linewidth=2.0)
end
plot([2.0],[2.0],"sk",linewidth=6.0)
plot([1.0],[0.0],"ok",linewidth=6.0)
annotate("\$(2,2)^\\top\$",[2.1,2.1],fontsize=18)
annotate("\$(1,0)^\\top\$",[1.25,0.05],fontsize=18)

ax = gca()
ax[:spines]["left"][:set_position]("center")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
xticks([]);yticks([])
xlim((-2xl,2xl));ylim((-2yl,2yl))
savefig(math2160*"/normAinfannotated.pdf";bbox_inches="tight")



clf();axes(aspect="equal")

p = 6
x = Float64[0.5 + i/2^p for i in 0:2^(p-1)-1]
z = zero(x)

for p in -6:3
    plot(2.0^p*x,z,"+k",markersize=p+7,markeredgewidth=2.0^(p-3))
end

ax = gca()
ax[:spines]["left"][:set_color]("none")
ax[:spines]["right"][:set_color]("none")
ax[:spines]["bottom"][:set_position]("center")
ax[:spines]["top"][:set_color]("none")
ax[:xaxis][:set_ticks_position]("bottom")
xticks([0;1/8;1/4;1/2;1;2;4],("\$0\$","\$\\frac{1}{8}\$","\$\\frac{1}{4}\$","\$\\frac{1}{2}\$","\$1\$","\$2\$","\$4\$"));yticks([])
xlim((0,6));ylim((-0.2,0.2))
savefig(math2160*"/FloatingPoint.pdf";bbox_inches="tight")


f = x -> log(1+x)/x
f1 = x -> log(1+x)/((1+x)-1)
f2 = x -> log1p(x)/x

xx = logspace(-16,0.0,340)
xxBF = logspace(big(-16),big(0.0),340)
fs = Float64[f(x) for x in xx];
fBFs = Float64[f(x) for x in xxBF];
f1s = Float64[f1(x) for x in xx];
f2s = Float64[f2(x) for x in xx];

err = abs.((fs-fBFs)./fBFs)
err1 = abs.((f1s-fBFs)./fBFs)
err2 = abs.((f2s-fBFs)./fBFs)

clf();
loglog(xx,err,"-r",xx,err1,"-g",xx,err2,"-b")
xlabel("\$x\$");ylabel("Relative Error");grid(true)
ylim((1e-17,1))
legend(["Naïve","First Alternative","Second Alternative"],loc="upper right")
savefig(math2160*"/log1pdx.pdf")


function sumjk(j,k)
    ret = 0.0
    for i=j:k
        ret += (-1.0)^i/i
    end
    ret
end

function sumjkrev(j,k)
    ret = 0.0
    for i=k:-1:j
        ret += (-1.0)^i/i
    end
    ret
end

function sumjksplit(j,k)
    t1 = 0.0
    for i=j:2:k
        t1 += inv(i)
    end
    t2 = 0.0
    for i=j+1:2:k
        t2 += inv(i)
    end
    if isodd(j)
        return t2-t1
    else
        return t1-t2
    end
end

function sumjkbf(j,k)
    ret = big(0.0)
    for i=big(j):big(k)
        ret += (big(-1.0))^i/i
    end
    ret
end

function sumpairwise(x::Vector)
    if length(x) == 1
        return x[1]
    else
        n2 = floor(Int,length(x)/2)
        return sumpairwise(x[1:n2])+sumpairwise(x[n2+1:end])
    end
end
sumjkpairwise(j,k) = sumpairwise((-1.0).^(j:k)./(j:k))

kk = round.([Int], logspace(2,6,170))

Sjk = Float64[sumjk(1,k) for k in kk]
Sjkrev = Float64[sumjkrev(1,k) for k in kk]
Sjksplit = Float64[sumjksplit(1,k) for k in kk]
Sjkpairwise = Float64[sumjkpairwise(1,k) for k in kk]
Sjkbf = Float64[sumjkbf(1,k) for k in kk]

errfor = abs.((Sjk-Sjkbf)./Sjkbf)
errrev = abs.((Sjkrev-Sjkbf)./Sjkbf)
errsplit = abs.((Sjksplit-Sjkbf)./Sjkbf)
errpairwise = abs.((Sjkpairwise-Sjkbf)./Sjkbf)

clf();
loglog(kk,errfor,"-r",kk,errrev,"-g",kk,errpairwise,"-b")
xlabel("\$n\$");ylabel("Relative Error");grid(true)
ylim((1e-17,1e-11))
legend(["Forward Summation","Reverse Summation","Pairwise Summation"],loc="upper right")
savefig(math2160*"/summation.pdf")

function bisection(a, b, f, tol)
    fa = f(a)
    fb = f(b)
    if fa*fb > 0
        error("Bisection cannot guarantee a root")
    end
    if fa == 0
        return a
    elseif fb == 0
        return b
    end
    absbma = abs(b-a)
    while abs(b-a) > tol*absbma
        c = (a+b)/2
        fc = f(c)
        if fa*fc > 0
           a = c
           fa = fc
        elseif fa*fc < 0
           b = c
           fb = fc
        else
           return c
        end
    end
    c = (a+b)/2
    return c
end

function bisection2(a, b, f, tol)
    fa = f(a)
    fb = f(b)
    if fa*fb > 0
        error("Bisection cannot guarantee a root")
    end
    if fa == 0
        return a
    elseif fb == 0
        return b
    end
    nmax = ceil(Int, log2(abs(b-a)/tol))
    c = (a+b)/2
    for n = 1:nmax
        c = (a+b)/2
        fc = f(c)
        if fa*fc > 0
           a = c
           fa = fc
        elseif fa*fc < 0
           b = c
           fb = fc
        else
           return c
        end
    end
    return c
end

f = x -> 2*cos(2x) + x^3 - exp(x)
c1 = bisection(0, 5, f, eps())
c2 = bisection2(0, 5, f, eps())

g = x -> sin(15x)
c1 = bisection(-1, 1, g, eps())
c2 = bisection2(-1, 1, g, eps())


hilbert(n::Int) = 1./((1:n).+(1:n)'-1)

function hilbertinv(n::Int)
    Hi = zeros(n,n)
    Hi[1] = n^2
    for i=2:n
        Hi[i,i] = (2i-3)*((n-i+1)*(n+i-1))^2//((i-1)^4*(2i-1))*Hi[i-1,i-1]
    end
    for i=2:n,j=1:i-1
        Hi[i,j] = -(i+j-2)*(n-i+1)*(n+i-1)//((i-1)^2*(i+j-1))*Hi[i-1,j]
    end
    for j=2:n,i=1:j-1
        Hi[i,j] = -(i+j-2)*(n-j+1)*(n+j-1)//((j-1)^2*(i+j-1))*Hi[i,j-1]
    end
    Hi
end

function hilbertinvtest(n::Int)
    Hi = zeros(BigFloat,n,n)
    for i=1:n,j=1:n
        Hi[i,j] = (-1.0)^(i+j)*(i+j-1)*binomial(n+i-big(1),n-j)*binomial(n+j-big(1),n-i)*binomial(i+j-big(2),i-1)^2
    end
    Hi
end

cond2Hn = [norm(hilbert(n),2)*norm(hilbertinv(n),2) for n in 5:5:200]

szego(n::Int) = (sqrt(2)+1)^(4n+4)/(2^(15/4)*sqrt(π*n))
cond2Hnasy = [szego(n) for n in 5:5:200]
@show [5:5:200 cond2Hn cond2Hnasy]
semilogy(5:5:200,cond2Hn,"-r",5:5:200,cond2Hnasy,"-g")
