import Pkg

Pkg.add(["FFTW", "PyPlot", "SpecialFunctions"])

using FFTW, LinearAlgebra, PyPlot, SparseArrays, SpecialFunctions

dftmatrix = N -> exp.(-2π*im*(0:N-1)*(0:N-1)'/N)

function permutationmatrix(N)
    p = [1:2:N-1;2:2:N]
    Matrix{Float64}(I, N, N)[:,p]
end

function Ω(N)
    ω = exp(-π*im/N)
    diagm(ω.^(0:N-1))
end

function negΩ(N)
    ω = exp(-π*im/N)
    diagm(ω.^(N+(0:N-1)))
end

function testfactorization(N)
    IN2 = Matrix{Float64}(I, N÷2, N÷2)
    FN = dftmatrix(N)
    FN2 = dftmatrix(N÷2)
    PN = permutationmatrix(N)
    norm(FN*PN - [IN2 Ω(N÷2); IN2 negΩ(N÷2)]*[FN2 zero(FN2); zero(FN2) FN2])
end

f = z -> exp(z)

N = 2^4
z = exp.(2π*im*(0:N-1)/N)
@show [real(fft(f.(z))/N) inv.(gamma.(1:N))]

dctmatrix = N -> 2cospi.((0:N-1)*(0.5:N-0.5)'/N)

D = dctmatrix(N)
norm(D*D'-diagm([4N;2N*ones(N-1)]))

N = 2^4
x = cospi.((0.5:N-0.5)/N)
fhatT = FFTW.r2r(f.(x),FFTW.REDFT10)/N; fhatT[1] *= 0.5;
ftrueT = besseli.(0:N-1,1); ftrueT[2:end] *= 2;
@show [fhatT ftrueT]

dstmatrix = N -> 2sinpi.((1:N)*(1:N)'/(N+1))

D = dstmatrix(N)
norm(D*D'-2(N+1)*I)

N = 2^4
x = cospi.((1:N)/(N+1))
wr = sinpi.((1:N)/(N+1))
fhatU = FFTW.r2r(wr.*f.(x), FFTW.RODFT00)/(N+1)
ftrueU = [1.1303182079849698,0.5429906790681531,0.13301054954599148,0.021896961768374933,0.0027146315595697186,0.0002698639377257711,2.2389055236813955e-5,1.5936998453382377e-6,9.933094552965612e-8,5.505896079673745e-9,2.747752278683482e-10,1.2469826768142837e-11,5.188642363338741e-13,1.9932612151559182e-14,7.111389153842272e-16,2.3682880915332792e-17]
@show [fhatU ftrueU]


function conversion(N)
    A = spzeros(N,N)
    A[1] = 1
    for i=2:N
        A[i,i] = 0.5
    end
    for i = 1:N-2
        A[i,i+2] -= 0.5
    end
    A
end

@show norm(fhatU-conversion(N)*fhatT, Inf)

@show norm(fhatT-conversion(N)\fhatU, Inf)
