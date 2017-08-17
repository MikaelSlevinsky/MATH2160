immutable HyperCube{T,d}
    left::Vector{T}
    right::Vector{T}
    function HyperCube{T,d}(left::Vector{T}, right::Vector{T}) where {T,d}
        @assert length(left) == length(right) == d
        new{T,d}(left, right)
    end
end

HyperCube{T}(left::Vector{T}, right::Vector{T}) = HyperCube{T, length(left)}(left, right)
HyperCube{T}(::Type{T}, d::Int) = HyperCube(-ones(T, d), ones(T, d))
HyperCube(d::Int) = HyperCube(Float64, d)

center(hypercube::HyperCube) = 0.5*(hypercube.left + hypercube.right)
Base.rand{T,d}(hypercube::HyperCube{T,d}) = (hypercube.right - hypercube.left).*rand(T,d) + hypercube.left

immutable Problem{d}
    f::Function
    hypercube::HyperCube{Float64,d}
end

immutable Solution{d}
    problem::Problem{d}
    allparents::Vector{Vector{Vector{Float64}}}
    allchildren::Vector{Vector{Vector{Float64}}}
    location::Vector{Float64}
    minimum::Float64
end

"""
    Implement an evolutionary search based on prokayotic (mitotic) reproduction.

    Input: a Julia object of the type Problem.

    Output: a Julia object of the type Solution.
"""
function prokaryoticsolve(problem::Problem; nmembers::Int = 50, ngenerations::Int = 25, withprint::Bool = true, scale::Float64 = 0.5)
    f, hypercube = problem.f, problem.hypercube

    generation, idxbest, h = 1, 1, 1.0

    parents = [center(hypercube) for i in 1:nmembers]
    fvals = map(f, parents)

    allparents = [copy(parents)]
    allchildren = Vector{Vector{Float64}}[]
    while generation ≤ ngenerations
        children = copy(parents)
        for p in parents
            append!(children, [p + h*(rand(hypercube) - center(hypercube)) for i in 2:nmembers])
        end
        fvals = map(f,children)
        idx = sortperm(fvals)
        copy!(parents, children[idx[1:nmembers]])
        push!(allparents, copy(parents))
        push!(allchildren, copy(children))
        idxbest = idx[1]
        h *= scale
        withprint && println("Generation: $(lpad(generation,ceil(Int,log10(ngenerations)))), and the fitness of the best survivor: ", fvals[idxbest])
        generation += 1
    end
    Solution(problem, allparents, allchildren, parents[1], fvals[idxbest])
end

"""
    Implement an evolutionary search based on eukayotic (meiotic) reproduction.

    Input: a Julia object of the type Problem.

    Output: a Julia object of the type Solution.
"""
function eukaryoticsolve(problem::Problem; nmembers::Int = 50, ngenerations::Int = 25, withprint::Bool = true)
    f, hypercube = problem.f, problem.hypercube

    generation, idxbest = 1, 1

    parents = [rand(hypercube) for i in 1:nmembers]
    fvals = map(f, parents)

    allparents = [copy(parents)]
    allchildren = Vector{Vector{Float64}}[]

    while generation ≤ ngenerations
        children = Vector{Float64}[]
        for p in parents
            append!(children, [p + rand()*(q-p) for q in parents])
        end
        fvals = map(f,children)
        idx = sortperm(fvals)
        copy!(parents, children[idx[1:nmembers]])
        push!(allparents, copy(parents))
        push!(allchildren, copy(children))
        idxbest = idx[1]
        withprint && println("Generation: $(lpad(generation,ceil(Int,log10(ngenerations)))), and the fitness of the best survivor: ", fvals[idxbest])
        generation += 1
    end
    Solution(problem, allparents, allchildren, parents[1], fvals[idxbest])
end


################################################################################
# Sample 2D function from course notes.
################################################################################

f = x -> exp(sin(50x[1])) + sin(60exp(x[2])) + sin(70sin(x[1])) + sin(sin(80x[2])) - sin(10*(x[1]+x[2])) + (x[1].^2 + x[2].^2)./4

problem = Problem(f, HyperCube(2))

solution = prokaryoticsolve(problem);

solution = eukaryoticsolve(problem);

################################################################################
# The Rosenbrock function that is considered hard to minimize.
################################################################################

f = x -> (1-x[1])^2 + 100*(x[2]-x[1]^2)^2

problem = Problem(f, HyperCube([-2.,-2.],[2.,2.]))

solution = prokaryoticsolve(problem);

solution = eukaryoticsolve(problem);
