#!/usr/bin/env julia
# gMAP3bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to this
# time done for the reduced 3 species case

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
using SymEngine # Trying alternative substitution method to see if this is faster
import GR # Need this to stop world age plotting error?

# first should make generic functions that take in a full set of parameters, that
# can perform the machine alegbra and then sub in the parameter set to generate
# numeric values for each of the instances considered

# make a symbolic diffusion matrix
function Ds()
    A, B, S, K, k, Q, q, kmin, qmin, f, r, F = symbols("A B S K k Q q kmin qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(3,3)
    e[1,2:3] = e[2,1] = e[2,3] = 0 # maybe e[3,3] here
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    A, B, S, K, k, Q, q, kmin, qmin, f, r, F = symbols("A B S K k Q q kmin qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(3)
    b[1] = k*S*r/(r + f*B^2) - kmin*A - K*A
    b[2] = q*S*r/(r + f*A^2) - qmin*B - Q*B
    b[3] = -k*S*r/(r + f*B^2) - q*S*r/(r + f*A^2) + kmin*A + qmin*B + F
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    the1, the2, the3 = symbols("the1 the2 the3")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the1*b[1] + the2*b[2] + the3*b[3]
    H += 0.5*the1*(D[1,1]*the1 + D[1,2]*the2 + D[1,3]*the3)
    H += 0.5*the2*(D[2,1]*the1 + D[2,2]*the2 + D[2,3]*the3)
    H += 0.5*the3*(D[3,1]*the1 + D[3,2]*the2 + D[3,3]*the3)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    A, B, S = symbols("A B S")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{SymEngine.Basic,1}(3)
    Hx[1] = diff(H, A)
    Hx[2] = diff(H, B)
    Hx[3] = diff(H, S)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    the1, the2, the3 = symbols("the1 the2 the3")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{SymEngine.Basic,1}(3)
    Hθ[1] = diff(H, the1)
    Hθ[2] = diff(H, the2)
    Hθ[3] = diff(H, the3)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    the1, the2, the3 = symbols("the1 the2 the3")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{SymEngine.Basic,2}(3,3)
    Hθθ[1,1] = diff(Hθ[1],the1)
    Hθθ[1,2] = diff(Hθ[1],the2)
    Hθθ[1,3] = diff(Hθ[1],the3)
    Hθθ[2,1] = diff(Hθ[2],the1)
    Hθθ[2,2] = diff(Hθ[2],the2)
    Hθθ[2,3] = diff(Hθ[2],the3)
    Hθθ[3,1] = diff(Hθ[3],the1)
    Hθθ[3,2] = diff(Hθ[3],the2)
    Hθθ[3,3] = diff(Hθ[3],the3)
    return(Hθθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ follwed by X
function Hθxs()
    A, B, S = symbols("A B S")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{SymEngine.Basic,2}(3,3)
    Hθx[1,1] = diff(Hθ[1],A)
    Hθx[1,2] = diff(Hθ[1],B)
    Hθx[1,3] = diff(Hθ[1],S)
    Hθx[2,1] = diff(Hθ[2],A)
    Hθx[2,2] = diff(Hθ[2],B)
    Hθx[2,3] = diff(Hθ[2],S)
    Hθx[3,1] = diff(Hθ[3],A)
    Hθx[3,2] = diff(Hθ[3],B)
    Hθx[3,3] = diff(Hθ[3],S)
    return(Hθx)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y1, y2, y3 = symbols("y1 y2 y3")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1] + b[1]*Dmin[1,2]*b[2] + b[1]*Dmin[1,3]*b[3]
    num += b[2]*Dmin[2,1]*b[1] + b[2]*Dmin[2,2]*b[2] + b[2]*Dmin[2,3]*b[3]
    num += b[3]*Dmin[3,1]*b[1] + b[3]*Dmin[3,2]*b[2] + b[3]*Dmin[3,3]*b[3]
    den = 0
    den += y1*Dmin[1,1]*y1 + y1*Dmin[1,2]*y2 + y1*Dmin[1,3]*y3
    den += y2*Dmin[2,1]*y1 + y2*Dmin[2,2]*y2 + y2*Dmin[2,3]*y3
    den += y3*Dmin[3,1]*y1 + y3*Dmin[3,2]*y2 + y3*Dmin[3,3]*y3
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2, y3 = symbols("y1 y2 y3")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(3)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    c[3] = λ*y3 - b[3]
    ϑ = Array{SymEngine.Basic,1}(3)
    ϑ[1] = Dmin[1,1]*c[1] + Dmin[1,2]*c[2] + Dmin[1,3]*c[3]
    ϑ[2] = Dmin[2,1]*c[1] + Dmin[2,2]*c[2] + Dmin[2,3]*c[3]
    ϑ[3] = Dmin[3,1]*c[1] + Dmin[3,2]*c[2] + Dmin[3,3]*c[3]
    return(ϑ)
end

# function to generate one of each symbolic object and to substitute parameter values into it
# this is a computationally intensive function, need to minimize the number of calls
function gensyms(ps::AbstractVector)
    # create symbolic objects
    ϑ = ϑs()
    λ = λs()
    Hθ = Hθs()
    Hθx = Hθxs()
    Hθθ = Hθθs()
    Hx = Hxs()
    H = Hs()
    # specify symbols that will be substituted for
    K, k, Q, q, kmin, qmin, f, r, F = symbols("K k Q q kmin qmin f r F")
    # now perform substitutions
    λ = subs(λ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
    λ = subs(λ, qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
    H = subs(H, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
    H = subs(H, qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
    for i = 1:3
        ϑ[i] = subs(ϑ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        ϑ[i] = subs(ϑ[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
        Hθ[i] = subs(Hθ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        Hθ[i] = subs(Hθ[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
        Hx[i] = subs(Hx[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        Hx[i] = subs(Hx[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
    end
    for j = 1:3
        for i = 1:3
            Hθx[i,j] = subs(Hθx[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
            Hθx[i,j] = subs(Hθx[i,j], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
            Hθθ[i,j] = subs(Hθθ[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
            Hθθ[i,j] = subs(Hθθ[i,j], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9])
        end
    end
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to find the zeros of the function
function nullcline(F::Number,r::Number,f::Number,K::Number,Q::Number,k::Number,q::Number,kmin::Number,qmin::Number,high2low::Bool)
    g(x) = (K - kmin)*(q/k)*((r + f*((F - K*x)/Q)^2)/(r + f*x^2))*x - (qmin + Q)*(F - K*x)/Q
    three = false
    n = 0
    As = []
    while three == false
        As = fzeros(g, 0, 2*F/K, order = 1)
        n = length(As)
        if n == 3
            three = true
        end
    end
    Bs = zeros(n)
    Ss = zeros(n)
    for i = 1:n
        Bs[i] = (F - K*As[i])/Q
        Ss[i] = (1/(k*r))*(r + f*((F - K*As[i])/Q)^2)*(K - kmin)*As[i]
    end
    sad = [ As[2]; Bs[2]; Ss[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ss[1] ]
        ss2 = [ As[3]; Bs[3]; Ss[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ss[3] ]
        ss2 = [ As[1]; Bs[1]; Ss[1] ]
    end
    print("$(ss1)\n")
    print("$(sad)\n")
    print("$(ss2)\n")
    return (ss1,sad,ss2)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray,NG::Int,Nmid::Int)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s1 = zeros(Nmid)
    s2 = zeros(NG+2-Nmid)
    s1[1] = 0
    s2[1] = 0
    for i = 2:Nmid
        dA = x[i,1] - x[i-1,1]
        dB = x[i,2] - x[i-1,2]
        dS = x[i,3] - x[i-1,3]
        s1[i] = s1[i-1] + sqrt(dA^2 + dB^2 + dS^2) # Could probably drop the sqrts to speed up the code
    end
    for i = 2:(NG+2-Nmid)
        dA = x[i+Nmid-1,1] - x[i+Nmid-2,1]
        dB = x[i+Nmid-1,2] - x[i+Nmid-2,2]
        dS = x[i+Nmid-1,3] - x[i+Nmid-2,3]
        s2[i] = s2[i-1] + sqrt(dA^2 + dB^2 + dS^2) # Could probably drop the sqrts to speed up the code
    end
    # Divide total arc length into equal segments
    ls1 = zeros(Nmid)
    ls2 = zeros(NG+2-Nmid)
    for i = 1:Nmid
        ls1[i] = (i-1)*s1[end]/(Nmid-1)
    end
    for i = 1:(NG+2-Nmid)
        ls2[i] = (i-1)*s2[end]/(NG+1-Nmid)
    end
    # Find first index greater than a ls[i] for each i
    inds1 = fill(0,Nmid)
    j = 1
    for i = 1:Nmid
        higher = false
        while higher == false
            if s1[j] >= ls1[i] || j == Nmid
                inds1[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    inds2 = fill(0,NG+2-Nmid)
    j = 1
    for i = 1:(NG+2-Nmid)
        higher = false
        while higher == false
            if s2[j] >= ls2[i] || j == NG + 2 - Nmid
                inds2[i] = j + Nmid - 1
                higher = true
            else
                j += 1
            end
        end
    end
    # First do mid points and end points as they should be fixed
    disx = zeros(NG+1,3)
    disx[1,:] = x[1,:]
    disx[Nmid,:] = x[Nmid,:]
    disx[NG+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:Nmid-1
        one = inds1[i] - 1
        two = inds1[i]
        s₀ = s1[one]
        s₁ = s1[two]
        for j = 1:3
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls1[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end

    for i = Nmid+1:NG
        one = inds2[i+1-Nmid] - 1
        two = inds2[i+1-Nmid]
        s₀ = s2[one+1-Nmid]
        s₁ = s2[two+1-Nmid]
        for j = 1:3
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

# function to generate the variables needed for a given algoritm iteration
function genvars(x::AbstractArray, λ::SymEngine.Basic, ϑ::Array{SymEngine.Basic,1}, NG::Int, Nmid::Int)
    # define neccesary symbols
    A, B, S, y1, y2, y3 = symbols("A B S y1 y2 y3")
    # calculate velocities
    xprim = fill(NaN, NG+1, 3)
    for i = 2:NG
        for j = 1:3
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 3)
    ϑt = Array{SymEngine.Basic,1}(3)
    for j = 1:3
        for i = 2:NG
            ϑt[j] = subs(ϑ[j], A=>x[i,1], B=>x[i,2], S=>x[i,3])
            ϑs[i,j] = subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3]) |> float
        end
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# function to be solved by NLsolve
function g!(F::AbstractArray, x::AbstractArray, C::AbstractVector, K::AbstractArray, xi::AbstractArray, NG::Int, Nmid::Int)
    for j = 1:3
        # Start point
        F[1,j] = x[1,j] - xi[1,j]
        # Mid point
        F[Nmid,j] = x[Nmid,j] - xi[2,j]
        # end point
        F[NG+1,j] = x[NG+1,j] - xi[3,j]
    end

    # first path
    for i = 2:Nmid-1
        for j = 1:3
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # second path
    for i = Nmid+1:NG
        for j = 1:3
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    return(F)
end

# function to solve the system of linear equations
function linsys(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,
                Hx::Array{SymEngine.Basic,1},Hθ::Array{SymEngine.Basic,1},Hθθ::Array{SymEngine.Basic,2},
                Hθx::Array{SymEngine.Basic,2},Δτ::Number,NG::Int,Nmid::Int,H::SymEngine.Basic)
    # define relevant symbols
    A, B, S, the1, the2, the3 = symbols("A B S the1 the2 the3")
    # Make array to store fixed points
    xi = fill(NaN, 3, 3)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Array{SymEngine.Basic,1}(3) # temporary hamiltonian so master isn't changed
    Hxt = Array{SymEngine.Basic,1}(3)
    for i = 1:3
        Hθt[i] = subs(Hθ[i], the1=>0.0, the2=>0.0, the3=>0.0)
        Hxt[i] = subs(Hx[i], the1=>0.0, the2=>0.0, the3=>0.0)
    end
    Hθtt = Array{SymEngine.Basic,1}(3)
    Hθttt = Array{SymEngine.Basic,1}(3)
    Ht = subs(H, the1=>0.0, the2=>0.0, the3=>0.0)
    Ht = subs(Ht, A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3]) |> float
    # loop to define fixed points
    for i = 1:3
        Hxt[i] = subs(Hxt[i], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3]) |> float
        Hθttt[i] = subs(Hθt[i], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3]) |> float
        Hθtt[i] = subs(Hθt[i], A=>x[1,1], B=>x[1,2], S=>x[1,3]) |> float
        Hθt[i] = subs(Hθt[i], A=>x[end,1], B=>x[end,2], S=>x[end,3]) |> float
        # start point
        xi[1,i] = Δτ*(Hθtt[i]) + x[1,i]
        # midpoint
        xi[2,i] = -Δτ*(Hxt[i]*Ht) + x[Nmid,i]
        # End point
        xi[3,i] = Δτ*(Hθt[i]) + x[end,i]
    end
    Hxθt = Array{SymEngine.Basic,2}(3,3)
    # loop to calculate Hxθ
    for j = 1:3
        for i = 1:3
            Hxθt[i,j] = subs(Hθx[i,j], the1=>0.0, the2=>0.0, the3=>0.0)
            Hxθt[i,j] = subs(Hxθt[i,j], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3]) |> float
        end
    end
    # now transpose to get right form
    Hxθt = transpose(Hxθt) # possible source of error if ive misunderstood here
    # loop for additional midpoint terms
    for j = 1:3
        for i = 1:3
            xi[2,i] -= Δτ*(Hxθt[i,j]*Hθt[j])
        end
    end

    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 3)
    Hxt = Array{SymEngine.Basic,1}(3)
    Hθθt = Array{SymEngine.Basic,2}(3,3)
    Hθxt = Array{SymEngine.Basic,2}(3,3)
    # Put initial values for K in
    for j = 1:3
        for i = 2:NG
            K[i,j] = x[i,j]
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end
    for l = 1:3
        for i = 2:NG
            # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
            Hxt[l] = subs(Hx[l], A=>x[i,1], B=>x[i,2], S=>x[i,3])
            Hxt[l] = subs(Hxt[l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3]) |> float
            for m = 1:3
                Hθθt[m,l] = subs(Hθθ[m,l], A=>x[i,1], B=>x[i,2], S=>x[i,3])
                Hθθt[m,l] = subs(Hθθt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3]) |> float
                Hθxt[m,l] = subs(Hθx[m,l], A=>x[i,1], B=>x[i,2], S=>x[i,3])
                Hθxt[m,l] = subs(Hθxt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3]) |> float
                # Update K's with new contributions from Hamiltonians
                K[i,m] -= Δτ*λs[i]*(Hθxt[m,l]*xprim[i,l])
                K[i,m] += Δτ*(Hθθt[m,l]*Hxt[l])
            end
        end
    end
    # Make an initial guess of the path, as prior path
    newxi = x
    # make f! a closure of g! for specific xi, C, K
    f!(F,x) = g!(F,x,C,K,xi,NG,Nmid)
    # Then put all this into the solver
    newx = nlsolve(f!, newxi)
    return(newx)
end

# main function generates symbolic matrices that the 3 variables can be subbed into
function gMAP(K,k,Q,q,kmin,qmin,f,r,F,NM::Int,NG::Int,Nmid::Int,Δτ,high2low::Bool)
    # make vector of these optimisation parameters
    paras = [ K; k; Q; q; kmin; qmin; f; r; F ]
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(paras)

    # Now generate an initial path to optimize over
    ss1, sad, ss2 = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low)
    a1 = collect(linspace(ss1[1],sad[1],Nmid))
    a2 = collect(linspace(sad[1],ss2[1],NG+2-Nmid))
    a = vcat(a1,a2[2:length(a2)])
    b1 = collect(linspace(ss1[2],sad[2],Nmid))
    b2 = collect(linspace(sad[2],ss2[2],NG+2-Nmid))
    b = vcat(b1,b2[2:length(b2)])
    s1 = collect(linspace(ss1[3],sad[3],Nmid))
    s2 = collect(linspace(sad[3],ss2[3],NG+2-Nmid))
    s = vcat(s1,s2[2:length(s2)])
    x = hcat(a,b,s)

    # Then appropriatly discretise the path such that it works with this algorithm
    x = discretise(x,NG,Nmid)

    # Set up method to tell if is converged
    convrg = false
    l = 0
    gr()
    xold = x
    while convrg == false
        x, xprim, λs, ϑs, λprim = genvars(x,λ,ϑ,NG,Nmid)
        newx = linsys(x,xprim,λs,ϑs,λprim,Hx,Hθ,Hθθ,Hθx,Δτ,NG,Nmid,H)
        xn = discretise(newx.zero,NG,Nmid)
        # delta is the sum of the differences of all the points in the path
        δ = 0
        for i = 1:NG+1
            for j = 1:3
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        S = Ŝ(xn,xprim,λs,ϑs,λprim,NG)
        print("$(δ),$(sum(S))\n")
        if l % 50 == 0
            plot(x[:,1],x[:,2])
            savefig("../Results/GraphAB$(l).png")
            plot(x[:,3])
            savefig("../Results/GraphS$(l).png")
            plot(S)
            savefig("../Results/S$(l).png")
            plot(ϑs)
            savefig("../Results/vartheta$(l).png")
        end
        l += 1
        # Now overwrite old x
        x = xn
        if δ <=  0.00025
            convrg = true
            print("$(l) steps to converge\n")
        end
    end
    return(x)
end

# Function to calculate the action of a given path
function Ŝ(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int)
    S = zeros(NG-1)
    for j = 1:3
        S[1] += (3/(2*NG))*(xprim[2,j]*ϑs[2,j])
        # Not excluding the midpoint as the contribution is vanishing
        for i = 3:NG-1
            S[i-1] += (1/NG)*(xprim[i,j]*ϑs[i,j])
        end
        S[NG-1] += (3/(2*NG))*(xprim[NG,j]*ϑs[NG,j])
    end
    return(S)
end

# function to find the times of each point
function times(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::AbstractVector,x::AbstractArray,NM::Int)
    # Make discrete vector of time points
    t = zeros(NM+1)
    for i = 1:NM+1
        t[i] = (i-1)*ts[end]/NM
    end
    # Find index of first element greater than t[i] in ts
    inds = fill(0,NM+1)
    j = 1
    for i = 1:NM+1
        higher = false
        while higher == false
            if ts[j] >= t[i] || j == NG+1
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they are fixed
    path = zeros(NM+1,3)
    path[1,:] = x[1,:]
    path[NM+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        for j = 1:3
            x₀ = x[one,j]
            x₁ = x[two,j]
            path[i,j] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
        end
    end
    return(path)
end

function main()
    # General parameters
    K = 1
    k = 1
    Q = 1
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 100 # Promoter switching
    r = 10
    F = 10


    # Optimisation parameters
    NM = 300 # number of segments to discretise MAP onto
    NG = 300 # number of segments to optimize gMAP over
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.01 # I've made this choice arbitarily, too large and the algorithm breaks
    high2low = false # Set if starting from high state or low state

    # Now call simulation function with these parameters
    path = gMAP(K,k,Q,q,kmin,qmin,f,r,F,NM,NG,Nmid,Δτ,high2low)
    plot(path[:,1],path[:,2])
    savefig("../Results/Graph1.png")
    plot(path[:,3])
    savefig("../Results/Graph2.png")
end


@time main()
