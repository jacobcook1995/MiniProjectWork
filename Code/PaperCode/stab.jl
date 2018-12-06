#!/usr/bin/env julia
# stab.jl
# A script to read in the steady states and parameters from other scripts and then
# calculate miniumum action paths between the two
#
# Author: Jacob Cook
# Date: December 2018

using SymEngine
using NLsolve
using Plots

# make a symbolic diffusion matrix
function Ds()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,2,2)
    e[1,2] = e[2,1] = 0
    e[1,1] = sqrt(k*r/(r + f*B^2) + kmin*A + K*A + Kmin)
    e[2,2] = sqrt(q*r/(r + f*A^2) + qmin*B + Q*B + Qmin)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    for j = 1:2
        for i = 1:2
            Dmin[i,j] = expand(Dmin[i,j])
        end
    end
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,2)
    b[1] = k*r/(r + f*B^2) - kmin*A - K*A + Kmin
    b[2] = q*r/(r + f*A^2) - qmin*B - Q*B + Qmin
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    the1, the2 = symbols("the1 the2")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the1*b[1] + the2*b[2]
    H += 0.5*the1*(D[1,1]*the1 + D[1,2]*the2)
    H += 0.5*the2*(D[2,1]*the1 + D[2,2]*the2)
    H = expand(H)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    A, B = symbols("A B")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{SymEngine.Basic,1}(undef,2)
    Hx[1] = diff(H, A)
    Hx[2] = diff(H, B)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    the1, the2 = symbols("the1 the2")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{SymEngine.Basic,1}(undef,2)
    Hθ[1] = diff(H, the1)
    Hθ[2] = diff(H, the2)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    the1, the2 = symbols("the1 the2")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{SymEngine.Basic,2}(undef,2,2)
    Hθθ[1,1] = diff(Hθ[1],the1)
    Hθθ[1,2] = diff(Hθ[1],the2)
    Hθθ[2,1] = diff(Hθ[2],the1)
    Hθθ[2,2] = diff(Hθ[2],the2)
    return(Hθθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ follwed by X
function Hθxs()
    A, B = symbols("A B")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{SymEngine.Basic,2}(undef,2,2)
    Hθx[1,1] = diff(Hθ[1],A)
    Hθx[1,2] = diff(Hθ[1],B)
    Hθx[2,1] = diff(Hθ[2],A)
    Hθx[2,2] = diff(Hθ[2],B)
    return(Hθx)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y1, y2 = symbols("y1 y2")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1] + b[1]*Dmin[1,2]*b[2]
    num += b[2]*Dmin[2,1]*b[1] + b[2]*Dmin[2,2]*b[2]
    den = 0
    den += y1*Dmin[1,1]*y1 + y1*Dmin[1,2]*y2
    den += y2*Dmin[2,1]*y1 + y2*Dmin[2,2]*y2
    λ = sqrt(num/den)
    return(λ)
end

# function to find a symbolic expression for ϑ the conjugate momentum that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2 = symbols("y1 y2")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,2)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    ϑ = Array{SymEngine.Basic,1}(undef,2)
    ϑ[1] = Dmin[1,1]*c[1] + Dmin[1,2]*c[2]
    ϑ[2] = Dmin[2,1]*c[1] + Dmin[2,2]*c[2]
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
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("K k Q q kmin Kmin qmin Qmin f r")
    # now perform substitutions
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    λ = subs(λ, k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
    λ = subs(λ, Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    H = subs(H, k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
    H = subs(H, Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    for i = 1:2
        ϑ[i] = subs(ϑ[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        ϑ[i] = subs(ϑ[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        Hθ[i] = subs(Hθ[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        Hθ[i] = subs(Hθ[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        Hx[i] = subs(Hx[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        Hx[i] = subs(Hx[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    end
    for j = 1:2
        for i = 1:2
            Hθx[i,j] = subs(Hθx[i,j], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
            Hθx[i,j] = subs(Hθx[i,j], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
            Hθθ[i,j] = subs(Hθθ[i,j], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
            Hθθ[i,j] = subs(Hθθ[i,j], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        end
    end
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to generate the variables needed for a given algoritm iteration
function genvars(x::AbstractArray,λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int,Nmid::Int)
    # define neccesary symbols
    A, B, y1, y2 = symbols("A B y1 y2")
    # calculate velocities
    xprim = fill(NaN, NG+1, 2)
    for i = 2:NG
        for j = 1:2
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 2)
    ϑt = Array{SymEngine.Basic,1}(undef,2)
    for j = 1:2
        for i = 2:NG
            ϑt[j] = subs(ϑ[j], A=>x[i,1], B=>x[i,2])
            ϑs[i,j] = subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2]) |> float
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
function g!(F::AbstractArray,x::AbstractArray,C::AbstractVector,K::AbstractArray,xi::AbstractArray,NG::Int64,Nmid::Int64)
    # Start point
    F[1,1] = x[1,1] - xi[1,1]
    F[1,2] = x[1,2] - xi[1,2]
    # first path
    for i = 2:Nmid-1
        for j = 1:2
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # midpoint
    F[Nmid,1] = x[Nmid,1] - xi[2,1]
    F[Nmid,2] = x[Nmid,2] - xi[2,2]
    # second path
    for i = Nmid+1:NG
        for j = 1:2
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # end point
    F[NG+1,1] = x[NG+1,1] - xi[3,1]
    F[NG+1,2] = x[NG+1,2] - xi[3,2]
    return(F)
end

# function to solve the system of linear equations
function linsys(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},
                Hx::Array{SymEngine.Basic,1},Hθ::Array{SymEngine.Basic,1},Hθθ::Array{SymEngine.Basic,2},
                Hθx::Array{SymEngine.Basic,2},Δτ::Float64,NG::Int64,Nmid::Int64,H::SymEngine.Basic)
    # define relevant symbols
    A, B, the1, the2 = symbols("A B the1 the2")
    # Make array to store fixed points
    xi = fill(NaN, 3, 2)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Array{SymEngine.Basic,1}(undef,2) # temporary hamiltonian so master isn't changed
    Hxt = Array{SymEngine.Basic,1}(undef,2)
    for i = 1:2
        Hθt[i] = subs(Hθ[i], the1=>0.0, the2=>0.0)
        Hxt[i] = subs(Hx[i], the1=>0.0, the2=>0.0)
    end
    Hθtt = Array{SymEngine.Basic,1}(undef,2)
    Hθttt = Array{SymEngine.Basic,1}(undef,2)
    Ht = subs(H, the1=>0.0, the2=>0.0)
    Ht = subs(Ht, A=>x[Nmid,1], B=>x[Nmid,2]) |> float
    # loop to define fixed points
    for i = 1:2
        Hxt[i] = subs(Hxt[i], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        Hθttt[i] = subs(Hθt[i], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        Hθtt[i] = subs(Hθt[i], A=>x[1,1], B=>x[1,2]) |> float
        Hθt[i] = subs(Hθt[i], A=>x[end,1], B=>x[end,2]) |> float
        # start point
        xi[1,i] = Δτ*(Hθtt[i]) + x[1,i]
        # midpoint
        xi[2,i] = -Δτ*(Hxt[i]*Ht) + x[Nmid,i] # x[Nmid,i] is changing between loops for some reason
        # End point
        xi[3,i] = Δτ*(Hθt[i]) + x[end,i]
    end
    Hxθt = Array{SymEngine.Basic,2}(undef,2,2)
    # loop to calculate Hxθ
    for j = 1:2
        for i = 1:2
            Hxθt[i,j] = subs(Hθx[i,j], the1=>0.0, the2=>0.0)
            Hxθt[i,j] = subs(Hxθt[i,j], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        end
    end
    # now transpose to get right form
    Hxθt = transpose(Hxθt) # possible source of error if ive misunderstood here
    # loop for additional midpoint terms
    for j = 1:2
        for i = 1:2
            xi[2,i] -= Δτ*(Hxθt[i,j]*Hθttt[j]) # maybe Hθttt[j]?
        end
    end
    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 2)
    Hxt = Array{SymEngine.Basic,1}(undef,2)
    Hθθt = Array{SymEngine.Basic,2}(undef,2,2)
    Hθxt = Array{SymEngine.Basic,2}(undef,2,2)
    # Put initial values for K in
    for j = 1:2
        for i = 2:NG
            K[i,j] = x[i,j]
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end
    for l = 1:2
        for i = 2:NG
            # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
            Hxt[l] = subs(Hx[l], A=>x[i,1], B=>x[i,2])
            Hxt[l] = subs(Hxt[l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
            for m = 1:2
                Hθθt[m,l] = subs(Hθθ[m,l], A=>x[i,1], B=>x[i,2])
                Hθθt[m,l] = subs(Hθθt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
                Hθxt[m,l] = subs(Hθx[m,l], A=>x[i,1], B=>x[i,2])
                Hθxt[m,l] = subs(Hθxt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
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

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray,NG::Int64,Nmid::Int64)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s1 = zeros(Nmid)
    s2 = zeros(NG+2-Nmid)
    s1[1] = 0
    s2[1] = 0
    for i = 2:Nmid
        dA = x[i,1] - x[i-1,1]
        dB = x[i,2] - x[i-1,2]
        s1[i] = s1[i-1] + sqrt(dA^2 + dB^2) # Could probably drop the sqrts to speed up the code
    end
    for i = 2:(NG+2-Nmid)
        dA = x[i+Nmid-1,1] - x[i+Nmid-2,1]
        dB = x[i+Nmid-1,2] - x[i+Nmid-2,2]
        s2[i] = s2[i-1] + sqrt(dA^2 + dB^2) # Could probably drop the sqrts to speed up the code
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
    disx = zeros(NG+1,2)
    disx[1,:] = x[1,:]
    disx[Nmid,:] = x[Nmid,:]
    disx[NG+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:Nmid-1
        one = inds1[i] - 1
        two = inds1[i]
        s₀ = s1[one]
        s₁ = s1[two]
        for j = 1:2
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
        for j = 1:2
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

# Function to generate an initial path then run the alorithm until a MAP is obtained
# This path is then returned for other functions
function gMAP(ps::Array{Float64,1},NG::Int64,Nmid::Int64,Δτ::Float64,ss1::Array{Float64,1},sad::Array{Float64,1},ss2::Array{Float64,1})
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(ps)
    # First find the steady states and saddle point
    a1 = collect(range(ss1[1],stop=sad[1],length=Nmid))
    a2 = collect(range(sad[1],stop=ss2[1],length=Nmid))
    a = vcat(a1,a2[2:end])
    b1 = collect(range(ss1[2],stop=sad[2],length=Nmid))
    b2 = collect(range(sad[2],stop=ss2[2],length=Nmid))
    b = vcat(b1,b2[2:end])
    x = hcat(a,b)
    # Then appropriatly discretise the path such that it works with this algorithm
    x = discretise(x,NG,Nmid)
    # Set up method to tell if is converged
    convrg = false
    l = 0
    while convrg == false
        x, xprim, λs, ϑs, λprim = genvars(x,λ,ϑ,NG,Nmid)
        newx = linsys(x,xprim,λs,ϑs,λprim,Hx,Hθ,Hθθ,Hθx,Δτ,NG,Nmid,H)
        xn = discretise(newx.zero,NG,Nmid)
        # delta is the sum of the differences of all the points in the path
        δ = 0
        for i = 1:NG+1
            for j = 1:2
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        l += 1
        if l % 200 == 0
            println("$l,$δ")
        end
        # Now overwrite old x
        x = xn
        if δ <= 0.000000005#0.00000000005
            convrg = true
            print("$(l) steps to converge\n")
        end
    end
    return(x)
end

function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 10
    ps = zeros(l,w)
    open(infile, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,w+1)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            for i = 1:w
                ps[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of steady states to be read
    infile = "../Results/Fig3Data/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    l = countlines(infile)
    w = 6
    steads = zeros(l,w)
    open(infile, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,w+1)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            for i = 1:w
                steads[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # preallocate to store paths
    NG = 600
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.001
    path1 = zeros(NG+1,2)
    path2 = zeros(NG+1,2)
    # run the stability analysis for each of the hundred steady states
    for i = 1:l
        println("Run number: $(i)")
        flush(stdout)
        path1 = gMAP(ps[i,:],NG,Nmid,Δτ,[steads[i,1],steads[i,2]],[steads[i,3],steads[i,4]],[steads[i,5],steads[i,6]])
        path2 = gMAP(ps[i,:],NG,Nmid,Δτ,[steads[i,5],steads[i,6]],[steads[i,3],steads[i,4]],[steads[i,1],steads[i,2]])
        # define sensible names for the output files
        outfile1 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2B.csv"
        outfile2 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2A.csv"
        # open files for writing
        out_file1 = open(outfile1, "w")
        for j = 1:size(path1,1)
            line = "$(path1[j,1]),$(path1[j,2])\n"
            write(out_file1, line)
        end
        close(out_file1)
        out_file2 = open(outfile2, "w")
        for j = 1:size(path2,1)
            line = "$(path2[j,1]),$(path2[j,2])\n"
            write(out_file2, line)
        end
        close(out_file2)
    end
    return(nothing)
end

@time main()
