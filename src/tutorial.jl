#
# Every programmer's introduction to a new language starts with a classic.
#

function main()
    println("Hello, world!")
end

#
# When we run `main()`, it calls the function println with the string "Hello, world!"
# and it appears in the console.
#

main()

#
# Defining our own functions couldn't be easier.
#

f = x -> exp(sin(x))

f(1)

exp(sin(1))

#
# Sometimes, however, the function body is more complicated (or more clearly
# written with more lines of code). Here is another example.We need to know
# how to write a function with conditional statements.
#

function myabs(x)
    if x < 0
        -x
    elseif x > 0
        x
    else
        0*x
    end
end

myabs(4)

myabs(-3.0)

#
# Another important set of structures are loops. Loops allow us to repeat a
# task programmatically.
#
# There are `while` loops and `for` loops. `While` loops continue execution while their
# conditional statement evaluates to the Boolean `true`. On the other hand, `for` loops
# execute by starting the iterable object, continuing through all its elements, and
# terminating at the end.
#

function myjourney(x)
    println("And $(x) goes")
    while true
        println("on and...")
    end
end

# Try `myjourney("it")`. To kill it, use `Ctrl+c`.

for k = 1:100
    println("This is k: ",k)
end

function myshorterjourney(x,k)
    print("And $(x) goes")
    i = 1
    while i ≤ k
        print(" on and")
        i += 1
    end
    println(" oooon. Strangers...")
end

# Try `myshorterjourney("it", 30)`

#
# Julia has full access to LAPACK (Linear Algebra PACKage), and thus can perform
# linear algebra very efficiently. In Julia, vectors and matrices are `subtypes`
# of the same abstract supertype: `Array`. An array is just like how it sounds:
# it is an indexable container for instances of a certain type. Let's give it a go.
#

using LinearAlgebra

# We create a vector of length 10 of all zeros:
x = zeros(10)

# We can update specific entries of `x` using the square brackets:

x[1] = 4
x[5] = 3
x[7] = 2

# We create a 10x10 matrix of all ones with twos on the diagonal:
A = ones(10, 10) + I

# Matrix-vector multiplication is native in Julia.
# Furthermore, the following creates the new vector `b`:

b = A*x

# Since the entries of `A` and `x` were integers (converted to floating-point
# numbers), `b` also consists of integers, which we check as:

b - Int.(b)

# Alternatively, we may look at the binary bits of every entry in `b`:

bitstring.(b)

# That strains the eyes, so let's use color to help us distinguish sign, exponent,
# and mantissa.

function colorbitstring(x::Float64)
    s = bitstring(x)
    print("\"")
    printstyled(string(s[1]); color = :red)
    printstyled(s[2:12]; color = :green)
    printstyled(s[13:end]; color = :blue)
    println("\"")
end

colorbitstring.(b);

# As they are terminating binary representations, we have reasonable evidence
# they are exactly computed.
# With `b` in hand, we could try to solve for `x` in `A*x = b`.
# The function `\`, invokes a polyalgorithm to solve the linear system
# `A*x = b` for `x`.
# Warning: while `inv(A)` calculates the matrix inverse A^{-1}, it is worse
# to do `inv(A)*b` than `A\b` for many reasons that we will discuss.

xapprox = A\b
err = x - xapprox

# We can now use norms to check the backward error, that is, the error in `xapprox`:

norm(err)

# To check a function for any documentation, type `?` to go into help mode.
# Then you can type the name of the function and press enter to check the docs.

# As we learned in class, the most important vector norms are the p-norms.
# These are implemented with another argument in the function signature.
# The one-argument `norm` falls back onto the 2-norm.

norm(err, 1)

norm(err, 2)

norm(err, 4)

norm(err, Inf)

# If `y` is a larger vector, then we cannot expect the backward error of
# `yapprox` to be of the same magnitude as that of `xapprox`. Therefore,
# numerical analysts are most often interested in the relative backward error:

norm(err, 1)/norm(x, 1)

norm(err, 2)/norm(x, 2)

norm(err, 4)/norm(x, 4)

norm(err, Inf)/norm(x, Inf)

# The residual is the error of `xapprox` in the image of `A`:

res = A*xapprox-b

norm(res)

# The norm of the residual divided by the norm of the approximate solution `x`
# is bounded by the condition number times `ulp`:

cond(A)

cond(A,1),cond(A,2),cond(A,Inf)

norm(res)/norm(x) ≤ cond(A)*eps()
norm(res,1)/norm(x,1) ≤ cond(A,1)*eps()
norm(res,Inf)/norm(x,Inf) ≤ cond(A,Inf)*eps()

# Induced matrix norms are also computed by a similar function:

opnorm(A)

# whereas the component-wise matrix norms are computed by the same function:

norm(A) # Frobenius norm

# The matrix we described visually in class was:

A = [1.0 2.0; 0.0 2.0]

opnorm(A, 1)

opnorm(A, 2)

opnorm(A, Inf)

# It's easy to manipulate vectors and matrices. For example:

A += 1e-3rand(2, 2)

# Similarly, suppose `b` were slightly different:

b[3] = sum(b)
b[8] = -9

# We can also check the inner-product (dot product) of `x` and `b`:

x'b
dot(x, b)
x⋅b

# We can also concatenate vectors:

y = [x;b]
y == vcat(x, b)

# We can extend this idea of concatenation to construction:

A = [1 2 3; 4 5 6; 7 8 9]
b = [1; 2; 3]

x = A\b

# Unfortunately, this matrix is singular:

det(A) == 0 # Note, not exactly 0 on all platforms

# and we got an appropriate warning.

# Some built-in functions help us with floating-point arithmetic.

x = floatmin()

colorbitstring(x)

x/2

colorbitstring(x/2)

x/2^52

colorbitstring(x/2^52)

(x/2^52)/2

x = floatmax()

colorbitstring(x)

colorbitstring(2x)

colorbitstring(2*(2x))

2x

# eps() is the relative unit of least precision. eps(1.0) == eps(),
# but eps(1500.0) > eps()

eps()

eps(1500.0) ≠ 1500eps()

# `eps(x)`` returns the spacing between `x` and the next
# representable floating-point number.

nextfloat(1500.0) - 1500

## Complex Numbers

# In Julia, complex numbers are just types with two fields to store the real
# and imaginary parts. The imaginary unit is `im`.

z = 1 + 2im

z.re
z.im

abs(z)

abs2(z)

z*conj(z)

w = 3 + 4im

w*z

# Let's try vectors with complex entries:

w = [-0.5+3im;2+2im]
z = [1.0+2im;3+4im]

dot(w, z)

norm(z) == sqrt(abs(z'z))

norm(z, 1)

norm(z, 2)

norm(z, Inf)

# Functions of a complex variable are similar to complex numbers. They can be
# represented as the sum of two separate functions as the purely real and purely
# imaginary parts. The exponential function of a purely imaginary number is
# particularly nice:

f = y -> exp(im*y)

g = y -> cos(y) + im*sin(y)

f(1)

g(1)

# There is an absolutely beautiful formula in mathematics
# that relates the five most important constants:

ℯ^(im*π) + 1

# Note that this italicized ℯ is obtained by \euler<TAB>

# For the exponential of a complex number,
#
# ℯ^z = ℯ^{x+iy} = ℯ^x(cos(y) + i sin(y))
#

f = z -> exp(z)

g = z -> exp(real(z))*(cos(imag(z)) + im*sin(imag(z)))

f(2+3im)

g(2+3im)


## Defining our own matrix operations couldn't be easier:

function matvec(A::Matrix, x::Vector)
    m, n = size(A)
    n == length(x) || throw(DimensionMismatch(
    "second dimension of A, $n, does not match length of x, $(length(x))"))
    b = zeros(m)
    for j = 1:n
        xj = x[j]
        for i = 1:m
            b[i] = b[i] + A[i,j]*xj
        end
    end
    b
end

# They perform almost as fast as the base library!

A = rand(1000, 100)
x = rand(100)

@time A*x;
@time matvec(A, x);

# Wait, first @time includes compile time.

@time A*x;
@time matvec(A, x);

# Ahh, much better. The second @time only includes run time.

# Functions are also useful in creating matrices,
# such as the Hilbert matrix in the assignment:

function hilbert(n)
    H = zeros(n, n)
    for j = 1:n
        for i = 1:n
            H[i,j] = inv(i+j-1)
        end
    end
    H
end

hilbert(15)

## Types and multiple dispatch

#
# In Julia, we can create objects and define their behaviour by how they
# interact with eachother and with functions. Objects in Julia are instances
# of (mutable) structs.
#

mutable struct Point2D
    x::Float64
    y::Float64
end

mutable struct Point3D
    x::Float64
    y::Float64
    z::Float64
end

# Here, the `::` notation is a type-assertion: x,y(, and z) must be 64-bit floating-point numbers.

p = Point2D(1.0, 2.0)
q = Point2D(3.0, 4.0)
r = Point3D(1.5, 2.5, 3.5)
s = Point3D(-2.0, 3.0, 0.5)

# `p`, `q`, `r`, and `s` are now instances of the Julia structs `Point2D` and `Point3D`.
# To access the fields, that is the coordinates, we follow the instance by a dot
# and the field in question:

p.x

# Julia functions support multiple dispatch: this means that we can define new methods
# of a function which only get called when a specific signature is present.

dist(p::Point2D, q::Point2D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2)
dist(p::Point3D, q::Point3D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2+(p.z-q.z)^2)

dist(p, q)
dist(r, s)

dist(p, s)

# It looks like a method is missing from the function `dist`. Let's add it:

dist(p::Point2D, q::Point3D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2+q.z^2)

# Great! but what about:

dist(r, q)

# Methods can call other methods:

dist(p::Point3D, q::Point2D) = dist(q, p)

dist(r, q)

# Here is another example of a type for the Maclaurin series of a function.

struct Maclaurin{T}
    coefficients::Vector{T}
end

ncoefficients(M::Maclaurin) = length(M.coefficients)

# We could use Horner's rule to evaluate the polynomial.

function horner(c::Vector, x)
    if isempty(c)
        return zero(x)
    end

    ret = zero(promote_type(eltype(c), eltype(x)))
    for k in length(c):-1:1
        ret = x*ret + c[k]
    end

    ret
end

(M::Maclaurin)(x) = horner(M.coefficients, x)

# Next, we construct a Maclaurin series to e^x. The syntax [... for i = 0:16]
# creates a comprehension. It's another shorthand for creating an array.

f = Maclaurin([1/factorial(i) for i = 0:16])

f(1) - exp(1)

f(-1) - exp(-1)

#
# To get access to Julia in Terminal from any folder on a mac,
# first check that you have a `.bash_profile`:
#
# cd ~
# ls -all
#
# If you do, execute this line:
#
# echo 'export PATH="$(pwd):/Applications/Julia-1.5.app/Contents/Resources/julia/bin:$PATH"' >> .bash_profile
#
# If you don't, create a `.bash_profile` via `touch .bash_profile` and execute the above.
#
