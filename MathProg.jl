const math2160 = "/Users/Mikael/Documents/Courses/MATH2160"

using PyPlot

import PyPlot: semilogy, contourf, surf

include("EvolutionaryAlgorithms.jl")

function semilogy(solution::Solution)
    nmembers, ngenerations = length(first(solution.allparents)), length(solution.allparents)
    clf()
    for n = 1:ngenerations
        semilogy(fill(n,nmembers), map(solution.problem.f, solution.allparents[n]) - solution.minimum,"o")
    end
    xlim((0,ngenerations))
    ylim((1e-16,1e1))
    xlabel("Generation")
    ylabel("Error")
    grid(true);gcf()
end

function plotprob(pop,prob,A,c)
    clf()
    plot(pop,prob,"ok")
    plot(pop,A*c,"-k",linewidth=2.0)
    xlim(extrema(pop))
    ylim((0,1))
    xlabel("Population")
    ylabel("Probability of success after \$25\$ generations")
    grid(true);gcf()
end

function contourf(problem::Problem{2},N=301,L=50;cmap="seismic")
    f, hypercube = problem.f, problem.hypercube

    xl, xr = hypercube.left[1], hypercube.right[1]
    yl, yr = hypercube.left[2], hypercube.right[2]

    x = linspace(xl, xr, N)
    y = linspace(yl, yr, N)
    xx, yy = [x for y in y, x in x], [y for y in y, x in x]
    fplotvals = map(f,map(tuple,xx,yy))

    clf()
    contourf(x,y,fplotvals,L;cmap=cmap)
    xlim((xl,xr))
    ylim((yl,yr))
    xlabel("\$x_1\$")
    ylabel("\$x_2\$")
    grid(false);gcf()
end

function surf(problem::Problem{2},N=301;cmap="seismic")
    f, hypercube = problem.f, problem.hypercube

    xl, xr = hypercube.left[1], hypercube.right[1]
    yl, yr = hypercube.left[2], hypercube.right[2]

    x = linspace(xl, xr, N)
    y = linspace(yl, yr, N)
    xx, yy = [x for y in y, x in x], [y for y in y, x in x]
    fplotvals = map(f,map(tuple,xx,yy))

    clf()
    surf(x,y,fplotvals;rstride=3, cstride=3, cmap=cmap, linewidth = 0.0,alpha=0.7)
    xlim((xl,xr))
    ylim((yl,yr))
    xlabel("\$x_1\$")
    ylabel("\$x_2\$")
    grid(false);gcf()
end

#=
cart2polr(x, y) = hypot(x, y)
cart2polt(x, y) = atan2(y, x)
@vectorize_2arg Number cart2polr
@vectorize_2arg Number cart2polt

function polarplot(solution::Solution)
    nmembers, ngenerations = length(first(solution.allparents)), length(solution.allparents)
    plot(;legend=false, proj = :polar)#, xscale = :log10, yscale=:log10)
    for n = 1:ngenerations
        x = map(x->x[1]-solution.location[1], solution.allparents[n])
        y = map(y->y[2]-solution.location[2], solution.allparents[n])
        if all(log10(cart2polr(x, y)) .!= 0.0)
            scatter!(cart2polt(x, y), 16 + log10(cart2polr(x, y)); markersize=0.5*(ngenerations+1-n))
        end
    end
    #ylims!(0.000001,16)
    current()
end

function makegif(solution::Solution{2})
    f, hypercube = solution.problem.f, solution.problem.hypercube

    xl, xr = hypercube.left[1], hypercube.right[1]
    yl, yr = hypercube.left[2], hypercube.right[2]

    x = linspace(xl, xr, 101)
    y = linspace(yl, yr, 101)
    xx, yy = [x for y in y, x in x], [y for y in y, x in x]
    fplotvals = map(f,map(tuple,xx,yy))

    contourf(x,y,fplotvals)
    xlims!(xl, xr)
    ylims!(yl, yr)

    allparents, allchildren = solution.allparents, solution.allchildren

    nmembers, ngenerations = length(first(allparents)), length(allparents)

    xp = map(x->x[1], allparents[1])
    yp = map(y->y[2], allparents[1])

    scatter!(xp, yp; legend=false, markersize = 1.0, color =:black)

    @gif for i = 2:ngenerations

        contourf(x,y,fplotvals)
        xlims!(xl, xr)
        ylims!(yl, yr)

        xc = map(x->x[1], allchildren[i-1])
        yc = map(y->y[2], allchildren[i-1])
        scatter!(xc, yc; legend=false, markersize = 1.0, color = :red)
        xp = map(x->x[1], allparents[i])
        yp = map(y->y[2], allparents[i])
        scatter!(xp, yp; legend=false, markersize = 1.0, color = :black)
    end
    current()
end
=#

################################################################################
# 2D
################################################################################

f = x -> exp(sin(50x[1])) + sin(60exp(x[2])) + sin(70sin(x[1])) + sin(sin(80x[2])) - sin(10*(x[1]+x[2])) + (x[1].^2 + x[2].^2)./4

problem = Problem(f, HyperCube(2))

solution = prokaryoticsolve(problem)

semilogy(solution)
#savefig(math2160*"/EvolutionaryAlgorithm1error.pdf")

solution = eukaryoticsolve(problem)

semilogy(solution)
#savefig(math2160*"/EvolutionaryAlgorithm2error.pdf")


#=
# Probability of success script. Takes a whole night.
nprokaryoticsuccesses = Vector{Int64}(130)
neukaryoticsuccesses = Vector{Int64}(130)

for nmembers = 1:130
    nprokaryoticsuccesses[nmembers] = 0
    neukaryoticsuccesses[nmembers] = 0
    for i = 1:1000
        solution = prokaryoticsolve(problem; nmembers = nmembers, withprint = false)
        solution.minimum < -3.3068686 && (nprokaryoticsuccesses[nmembers] += 1)
        solution = eukaryoticsolve(problem; nmembers = nmembers, withprint = false)
        solution.minimum < -3.3068686 && (neukaryoticsuccesses[nmembers] += 1)
    end
    println("nmembers = ",nmembers,". Number of prokaryotic successes: ",nprokaryoticsuccesses[nmembers]," and number of eukaryotic successes: ",neukaryoticsuccesses[nmembers])
end
=#

# Results
nprokaryoticsuccesses = [0,0,1,12,38,67,92,145,182,244,218,286,339,387,422,463,475,551,548,605,626,656,674,736,730,767,770,823,802,855,878,864,894,911,912,930,933,935,913,940,954,956,966,967,970,971,982,982,977,984,983,983,987,984,992,996,993,999,994,994,998,997,997,997,999,999,999,1000,1000,999,998,999,998,998,1000,999,999,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,999,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000]
neukaryoticsuccesses = [0,0,0,0,0,0,4,0,6,20,22,35,49,62,89,104,133,168,176,186,211,239,267,296,304,326,359,411,418,432,462,462,493,519,526,560,588,606,648,655,656,700,695,722,732,722,745,764,795,813,813,809,829,863,852,867,852,895,880,894,901,912,907,900,922,925,934,955,940,944,947,938,960,952,963,967,961,975,976,973,970,978,977,983,983,980,984,983,983,991,984,987,990,991,992,995,994,997,994,992,995,997,997,996,995,999,997,995,996,1000,1000,998,994,1000,999,998,999,999,998,999,1000,998,999,999,999,999,1000,999,1000,1000]

pop = 1:length(nprokaryoticsuccesses)
Pprokaryotic = nprokaryoticsuccesses/1000
Peukaryotic = neukaryoticsuccesses/1000

xpts = 2(pop-mean(pop))./(length(pop)-1) # Map 1:130 to [-1,1].

#=
# How many Chebyshev coefficients should I use?
# I balanced oscillations against positivity of a CDF.
for n in 1:30
    A = cos((0:n)'.*acos(xpts))

    c_prokaryotic = A\Pprokaryotic
    c_eukaryotic = A\Peukaryotic
    ext_pro = extrema(A*c_prokaryotic)
    ext_eu = extrema(A*c_eukaryotic)
    println("n = ",n," and: ",abs(ext_pro[1]) + abs(ext_pro[2]-1) + abs(ext_eu[1]) + abs(ext_eu[2]-1))
end
=#

A = cos((0:15)'.*acos(xpts))

c_prokaryotic = A\Pprokaryotic
c_eukaryotic = A\Peukaryotic

plotprob(pop,Pprokaryotic,A,c_prokaryotic)
#savefig(math2160*"/EvolutionaryAlgorithm1prob.pdf")
plotprob(pop,Peukaryotic,A,c_eukaryotic)
#savefig(math2160*"/EvolutionaryAlgorithm2prob.pdf")


problem1 = Problem(f, HyperCube([-20.,-20.],[20.,20.]))

surf(problem1,501)
#savefig(math2160*"/MathProgf1.png";dpi=300)

problem2 = Problem(f, HyperCube(2))

contourf(problem2)
#savefig(math2160*"/MathProgf2.png";dpi=300)

################################################################################
# 2D Newton iteration.
################################################################################

f = x -> exp(sin(50x[1])) + sin(60exp(x[2])) + sin(70sin(x[1])) + sin(sin(80x[2])) - sin(10*(x[1]+x[2])) + (x[1].^2 + x[2].^2)./4

function g(x)
    ret = fill(-10cos(10*(x[1]+x[2])), 2)
    ret[1] += 50.*cos(50x[1]).*exp(sin(50x[1])) + 70cos(x[1]).*cos(70sin(x[1])) + x[1]/2
    ret[2] += 60exp(x[2]).*cos(60exp(x[2])) + 80cos(80x[2]).*cos(sin(80x[2])) + x[2]/2
    ret
end

function H(x)
    ret = fill(100sin(10(x[1]+x[2])), 2, 2)
    ret[1,1] += ( (50cos(50x[1])).^2 - 50^2*sin(50x[1]) ).*exp(sin(50x[1])) + 1/2
    ret[1,1] -= (70cos(x[1]))^2.*sin(70sin(x[1])) + 70sin(x[1]).*cos(70sin(x[1]))

    ret[2,2] += 60exp(x[2]).*cos(60exp(x[2])) - (60exp(x[2])).^2.*sin(60exp(x[2])) + 1/2
    ret[2,2] -= (80cos(80x[2])).^2*sin(sin(80x[2])) + 80^2*sin(80x[2]).*cos(sin(80x[2]))

    ret
end

for i = 10:10:180
    solution = prokaryoticsolve(problem; nmembers = i, withprint = false)
    println("i = ", i,", and the fitness is: ", solution.minimum)
end

for i = 10:10:180
    solution = eukaryoticsolve(problem; nmembers = i, withprint = false)
    println("i = ", i,", and the fitness is: ", solution.minimum)
end

setprecision(33333)

x = [big(-0.02);big(0.21)]

while norm(g(x)) > log(eps(x[1]))^2*eps(x[1])
    x += H(x)\-g(x)
    println("x = [parse(BigFloat,\"",x[1],"\");")
    println("parse(BigFloat,\"",x[2],"\")]")
end
