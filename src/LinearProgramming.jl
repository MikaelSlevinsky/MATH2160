"""
`LinearProgramming` is a module that solves linear programs based on the simplex method.

# Examples
```jldoctest
julia> include("LinearProgramming.jl")
Main.LinearProgramming

julia> using LinearAlgebra, Main.LinearProgramming

julia> A = [1 -1 2 -3 2; -1 2 0 3 0]
2×5 Array{Int64,2}:
  1  -1  2  -3  2
 -1   2  0   3  0

julia> b = [1; 2]
2-element Array{Int64,1}:
 1
 2

julia> c = [1; 3; 2; 2; 4]
5-element Array{Int64,1}:
 1
 3
 2
 2
 4

julia> LP = LinearProgram(A, b, c)
LinearProgram{Rational{Int64}}:
  1//1   3//1  2//1   2//1  4//1  0//1
  1//1  -1//1  2//1  -3//1  2//1  1//1
 -1//1   2//1  0//1   3//1  0//1  2//1

julia> x = simplex!(LP)
The minimum of 13//3 is attained at x = Rational{Int64}[0//1, 0//1, 3//2, 2//3, 0//1].
5-element Array{Rational{Int64},1}:
 0//1
 0//1
 3//2
 2//3
 0//1

julia> c'x + LP.objective
0//1

julia> norm(A*x - b)
0.0
```
"""
module LinearProgramming

    import Base: eltype, show, size

    using LinearAlgebra

    export LinearProgram, simplex, simplex!,
           simplex_phase1!, simplex_phase2!

    """
    A `LinearProgram` stores a linear programming problem in standard form:

        minimize    cᵀx,
        subject to  A x = b,
        and         x ≥ 0.
    """
    mutable struct LinearProgram{T}
        A::Matrix{T}
        b::Vector{T}
        c::Vector{T}
        objective::T
        idxB::Vector{Int}
        idxN::Vector{Int}
        function LinearProgram{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}, objective::T, idxB::Vector{Int}, idxN::Vector{Int}) where T
            m, n = size(A)
            @assert length(b) == m
            @assert length(c) == n
            @assert length(idxB) == m
            @assert length(idxN) == n-m
            new(A, b, c, objective, idxB, idxN)
        end
    end

    LinearProgram(A::Matrix{T}, b::Vector{T}, c::Vector{T}, objective::T, idxB::Vector{Int}, idxN::Vector{Int}) where T = LinearProgram{T}(A, b, c, objective, idxB, idxN)
    LinearProgram(A::Matrix{T}, b::Vector{T}, c::Vector{T}, idxB::Vector{Int}) where T = LinearProgram(A, b, c, zero(T), idxB, setdiff(1:length(c),idxB))
    LinearProgram(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where T = LinearProgram(A, b, c, collect(1:length(b)))
    LinearProgram(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where T <: Integer = LinearProgram(convert(Matrix{Rational{T}},A), convert(Vector{Rational{T}},b), convert(Vector{Rational{T}},c))
    LinearProgram(A::Matrix{S}, b::Vector{T}, c::Vector{U}) where {S,T,U} = (V = promote_type(S,T,U); LinearProgram(convert(Matrix{V},A), convert(Vector{V},b), convert(Vector{V},c)))

    size(LP::LinearProgram) = size(LP.A)
    size(LP::LinearProgram,k::Integer) = size(LP.A,k)

    eltype(LP::LinearProgram{T}) where T = T

    function show(io::IO, LP::LinearProgram)
        println(io, summary(LP), ":")
        Base.print_matrix(io, Any[[LP.c'; LP.A] [LP.objective; LP.b]])
    end

    # Flip signs of any row of A and b for which bᵢ < 0.
    function flipsigns!(A::Matrix, b::Vector)
        m, n = size(A)
        for i = 1:m
            if b[i] < 0
                for j = 1:n
                    A[i,j] = -A[i,j]
                end
                b[i] = -b[i]
            end
        end
        A, b
    end

    function phase1setup!(AN::Matrix{T}, b::Vector{T}) where T
        m, n = size(AN)
        flipsigns!(AN, b)
        AB = Matrix{T}(I, m, m)
        A = [AN AB]
        c = zeros(T,m+n)
        idxB = collect(1+n:m+n)
        for i in idxB
            c[i] += 1
        end
        A, copy(b), c, idxB
    end


    function toreducedform!(LP::LinearProgram{T}) where T
        A, b, c, objective, idxB, idxN = LP.A, LP.b, LP.c, LP.objective, LP.idxB, LP.idxN
        m, n = size(LP)
        ID = Matrix{T}(I, m, m)
        AB = view(A,1:m,idxB)
        AN = view(A,1:m,idxN)
        LUF = lu(AB)
        AN[:] = LUF\AN
        b[:] = LUF\b
        AB[:] = ID

        for j in idxN
            for i in 1:m
                c[j] -= c[idxB[i]]*A[i,j]
            end
        end
        for i in 1:m
            objective -= c[idxB[i]]*b[i]
            c[idxB[i]] = zero(T)
        end
        LP.objective = objective

        LP
    end

    function chooseq(LP::LinearProgram)
        m, n = size(LP)
        c, idxN = LP.c, LP.idxN

        qj = 1
        q = idxN[qj]
        cq = c[q]
        for j = 1:n-m
            if c[idxN[j]] < cq
                cq = c[idxN[j]]
                qj = j
            end
        end
        qj, idxN[qj]
    end

    function choosep(LP::LinearProgram,qj,q)
        m, n = size(LP)
        A, b, idxB = LP.A, LP.b, LP.idxB

        pi = 1
        while A[pi,q] ≤ 0
            pi += 1
        end

        bpiApiq = b[pi]/A[pi,q]
        for i = pi:m
            if b[i]/A[i,q] < bpiApiq && A[i,q] > 0
                bpiApiq = b[i]/A[i,q]
                pi = i
            end
        end
        pi, idxB[pi]
    end

    function swappq!(LP::LinearProgram,pi,p,qj,q)
        idxB, idxN = LP.idxB, LP.idxN
        idxB[pi],idxN[qj] = idxN[qj],idxB[pi]
        LP
    end

    function checkreducedcosts(LP::LinearProgram)
        c, idxN = LP.c, LP.idxN
        ret = true
        for j in idxN
            ret *= c[j] ≥ 0
        end

        ret
    end

    function simplex_phase1!(LP::LinearProgram)
        Anew, bnew, cnew, idxBnew = phase1setup!(LP.A, LP.b)
        LPnew = LinearProgram(Anew, bnew, cnew, idxBnew)
        simplex_phase2!(LPnew)
        if LPnew.objective < 0
            error("Phase 1 complete. No basic feasible solution exists.")
        end
        LP.idxB[:] = LPnew.idxB
        LP.idxN[:] = setdiff(1:length(LP.c),LP.idxB)
        LP
    end

    function simplex_phase2!(LP::LinearProgram)
        toreducedform!(LP)
        while !checkreducedcosts(LP)
            qj, q = chooseq(LP)
            pi, p = choosep(LP, qj, q)
            swappq!(LP, pi, p, qj, q)
            toreducedform!(LP)
        end
        LP
    end

    function simplex!(LP::LinearProgram)
        simplex_phase1!(LP)
        simplex_phase2!(LP)
        x = zeros(eltype(LP), size(LP, 2))
        x[LP.idxB] .= LP.b
        println("The minimum of ", -LP.objective, " is attained at x = ", x, ".")
        x
    end

    simplex(LP::LinearProgram) = simplex!(deepcopy(LP))
end
