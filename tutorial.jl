#
# Every programmer's introduction to a new language starts with a classic.
#

function main()
    println("Hello, world!")
end

main()

#
# When we run `main()`, it calls the function println with the string "Hello, world!"
# and it appears in the console.
#


#
# Next, we need to know how to write a function with conditional statements.
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

for k=1:100
    println("This is k: ",k)
end

function myshorterjourney(x,k)
    print("And $(x) goes")
    i=1
    while i≤k
        print(" on and")
        i+=1
    end
    println(" oooon. Strangers...")
end

#
# Julia has full access to LAPACK, and thus can perform linear algebra very efficiently.
# In Julia, vectors and matrices are `subtypes` of the same abstract supertype: `Array`.
# An array is just like how it sounds: it is an indexable container for instances of
# a certain type. Let's give it a go.
#

# We create a vector of length 10 of all zeros:
x = zeros(10)

# We can update specific entries of `x` using the square brackets:

x[1] = 4
x[5] = 3
x[7] = 2

# We create a 10x10 matrix of all ones with twos on the diagonal:
A = ones(10,10) + I

# then we add a small perturbation:

A += 1e-3rand(10,10)

# Matrix-vector multiplication is native in Julia. Furthermore, the following creates
# the new vector `b`:

b = A*x

# Similarly, suppose `b` were slightly different:

b[3] = sum(b)
b[8] = -9

# Then, the notation `\`, will invoke a polyalgorithm to solve the linear system
# `Ax = b` for `x`.

x = A\b

# Certain functions are available for Matrices and vectors.

res = A*x-b

norm(res)

# calculates the `2-norm` of the vector of residuals.
# ? is Julia help. With another argument:

norm(res,1),norm(res,2),norm(res,4),norm(res,Inf)

# `norm` calculates the `p-norm` of a vector.
# It also calculates the induced matrix norms:

norm(A)

norm(A,1),norm(A,2),norm(A,Inf)

# The norm of the residual divided by the norm of the approximate solution `x`
# is bounded by the condition number times `ulp`:

cond(A)

cond(A,1),cond(A,2),cond(A,Inf)

norm(res)/norm(x) ≤ cond(A)*eps()
norm(res,1)/norm(x,1) ≤ cond(A,1)*eps()
norm(res,Inf)/norm(x,Inf) ≤ cond(A,Inf)*eps()

# We can also check the dot-product of `x` and `b`:

x⋅b
dot(x,b)

# It's equivalent to `x'*b`, except that it returns a number rather than an Array.

x'*b

# We can also concatenate vectors:

y = [x;b]
y == vcat(x,b)

# We can extend this idea of concatenation to construction:

A = [1 2 3; 4 5 6; 7 8 9]
b = [1; 2; 3]

x = A\b

# Unfortunately, this matrix is singular:

det(A) == 0

# and we got an appropriate warning.

## Types and multiple dispatch

#
# In Julia, we can create objects and define their behaviour by how they
# interact with eachother and with functions. Objects in Julia are instances
# of `type`s.
#

type Point2D
    x::Float64
    y::Float64
end

type Point3D
    x::Float64
    y::Float64
    z::Float64
end

# Here, the `::` notation is a type-assertion: x,y(, and z) must be floating-point numbers.

p = Point2D(1.0,2.0)
q = Point2D(3.0,4.0)
r = Point3D(1.5,2.5,3.5)
s = Point3D(-2.0,3.0,0.5)

# Julia functions support `multiple dispatch`: this means that we can define new methods
# of a function which only get called when the specific argument signature is present.

dist(p::Point2D,q::Point2D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2)
dist(p::Point3D,q::Point3D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2+(p.z-q.z)^2)

dist(p,q)
dist(r,s)

dist(p,s)

# It looks like a method is missing from the function `dist`. Let's add it:

dist(p::Point2D,q::Point3D) = sqrt((p.x-q.x)^2+(p.y-q.y)^2+q.z^2)

# Great! but what about:

dist(r,q)

# Methods can call other methods:

dist(p::Point3D,q::Point2D) = dist(q,p)

dist(r,q)


## Complex Numbers

# In Julia, complex numbers are just `types` with two fields to store the real
# and imaginary parts. The imaginary unit is `im`.

w = 1+2im
z = 3+4im
w*z

abs(z)
z.re
z.im

conj(z)

# Let's try vectors with complex entries:

w = [-0.5+3im;2+2im]
z = [1.0+2im;3+4im]

dot(w,z)
norm(z) == sqrt(abs(dot(z,z)))

norm(z,1),norm(z,2),norm(z,Inf)


## Defining our own matrix operations:

function matvec(A::Matrix, x::Vector)
    m,n = size(A)
    if n != length(x)
        throw(DimensionMismatch("A and x are not of conformable sizes"))
    end
    b = zeros(m)
    for j=1:n
        for i=1:m
            b[i] = b[i] + A[i,j]*x[j]
        end
    end

    b
end

# They perform almost as fast as the base library! (In how many dynamic programming
# languages is this type of performance even possible?!?!)

A = rand(1000,100)
x = rand(100)

@time A*x;
@time matvec(A,x);

# Wait, first @time includes compile time.

@time A*x;
@time matvec(A,x);

# Ahh, much better.
