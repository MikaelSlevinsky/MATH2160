module LinearProgramming

    import Base: size, show

    export LinearProgram, simplex, simplex!,
           simplex_phase1!, simplex_phase2!

    """
        A `LinearProgram` stores a linear programming problem in standard form:

        minimize    cᵀx,
        subject to  A x = b,
        and         x ≥ 0.
    """
    type LinearProgram{T}
        A::Matrix{T}
        b::Vector{T}
        c::Vector{T}
        objective::T
        idxB::Vector{Int64}
        idxN::Vector{Int64}
        function LinearProgram{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}, objective::T, idxB::Vector{Int64}, idxN::Vector{Int64}) where T
            m, n = size(A)
            @assert length(b) == m
            @assert length(c) == n
            @assert length(idxB) == m
            @assert length(idxN) == n-m
            new(A, b, c, objective, idxB, idxN)
        end
    end

    LinearProgram{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}, objective::T, idxB::Vector{Int64}, idxN::Vector{Int64}) = LinearProgram{T}(A, b, c, objective, idxB, idxN)
    LinearProgram{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}, idxB::Vector{Int64}) = LinearProgram(A, b, c, zero(T), idxB, setdiff(1:length(c),idxB))
    LinearProgram{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}) = LinearProgram(A, b, c, collect(1:length(b)))
    LinearProgram{T<:Integer}(A::Matrix{T}, b::Vector{T}, c::Vector{T}) = LinearProgram(convert(Matrix{Rational{T}},A), convert(Vector{Rational{T}},b), convert(Vector{Rational{T}},c))
    LinearProgram{S,T,U}(A::Matrix{S}, b::Vector{T}, c::Vector{U}) = (V = promote_type(S,T,U); LinearProgram(convert(Matrix{V},A), convert(Vector{V},b), convert(Vector{V},c)))

    Base.size(LP::LinearProgram) = size(LP.A)
    Base.size(LP::LinearProgram,k::Integer) = size(LP.A,k)

    function Base.show(io::IO,LP::LinearProgram;header::Bool=true)
        header && println(io,summary(LP),":")
        cobjAb = Any[[LP.c.';LP.A] [LP.objective;LP.b]]
        Base.showarray(io,cobjAb,false;header=false)
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

    function phase1setup!{T}(AN::Matrix{T}, b::Vector{T})
        m, n = size(AN)
        flipsigns!(AN, b)
        AB = eye(T, m)
        A = [AN AB]
        c = zeros(T,m+n)
        idxB = collect(1+n:m+n)
        for i in idxB
            c[i] += 1
        end
        A, copy(b), c, idxB
    end


    function toreducedform!{T}(LP::LinearProgram{T})
        A, b, c, objective, idxB, idxN = LP.A, LP.b, LP.c, LP.objective, LP.idxB, LP.idxN
        m, n = size(LP)
        ID = eye(T, m)
        AB = view(A,1:m,idxB)
        AN = view(A,1:m,idxN)
        LUF = lufact(AB)
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
    end

    function simplex{T}(LP::LinearProgram{T})
        LPc = simplex!(deepcopy(LP))
        x = zeros(T,size(LP, 2))
        LP.idxB[:] = LPc.idxB
        LP.idxN[:] = LPc.idxN
        x[LP.idxB] = LPc.b
        println("The minimum of ",-LPc.objective," is attained at x = ",x,".")
        x
    end

end
