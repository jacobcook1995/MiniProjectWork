#!/usr/bin/env julia
# stabS.jl
# A script to read in the steady states and parameters from other scripts and then
# calculate miniumum action paths between the two
# this is now done for the two species Schlögl model
#
# Author: Jacob Cook
# Date: December 2018

using SymEngine
using NLsolve
using Plots

# make a symbolic diffusion matrix
function Ds()
    X, k1, K1, k2, K2 = symbols("X k1 K1 k2 K2")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,1,1)
    e[1,1] = sqrt(k2*X^3 + K2*X^2 + K1*X + k1)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    Dmin[1,1] = expand(Dmin[1,1])
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    X, k1, K1, k2, K2 = symbols("X k1 K1 k2 K2")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,1)
    b[1] = -k2*X^3 + K2*X^2 - K1*X + k1
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    the = symbols("the")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the*b[1]
    H += 0.5*the*D[1,1]*the
    H = expand(H)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    X = symbols("X")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{SymEngine.Basic,1}(undef,1)
    Hx[1] = diff(H, X)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    the = symbols("the")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{SymEngine.Basic,1}(undef,1)
    Hθ[1] = diff(H, the)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    the = symbols("the")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{SymEngine.Basic,2}(undef,1,1)
    Hθθ[1,1] = diff(Hθ[1],the)
    return(Hθθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ follwed by X
function Hθxs()
    X = symbols("X")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{SymEngine.Basic,2}(undef,1,1)
    Hθx[1,1] = diff(Hθ[1],X)
    return(Hθx)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y = symbols("y")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1]
    den = 0
    den += y*Dmin[1,1]*y
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

# function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y = symbols("y")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,1)
    c[1] = λ*y - b[1]
    ϑ = Array{SymEngine.Basic,1}(undef,1)
    ϑ[1] = Dmin[1,1]*c[1]
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
    k1, K1, k2, K2 = symbols("k1 K1 k2 K2")
    # now perform substitutions
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    λ = subs(λ, k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    H = subs(H, k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    ϑ[1] = subs(ϑ[1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    Hθ[1] = subs(Hθ[1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    Hx[1] = subs(Hx[1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    Hθx[1,1] = subs(Hθx[1,1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    Hθθ[1,1] = subs(Hθθ[1,1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to generate the variables needed for a given path
function genvars(x::Array{Float64,1},λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int64,Nmid::Int64)
    # define neccesary symbols
    X, y = symbols("X y")
    # calculate velocities
    xprim = fill(NaN, NG+1)
    for i = 2:NG
        xprim[i] = (x[i+1] - x[i-1])/(2/NG)
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λs[i] = subs(λt, X=>x[i], y=>xprim[i]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λs[i] = subs(λt, X=>x[i], y=>xprim[i]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1)
    ϑt = Array{SymEngine.Basic,1}(undef,1)
    for i = 2:NG
        ϑt[1] = subs(ϑ[1], X=>x[i])
        ϑs[i] = subs(ϑt[1], y=>xprim[i]) |> float
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# function to be solved by NLsolve
function g!(F::Array{Float64,1},x::Array{Float64,1},C::Array{Float64,1},K::Array{Float64,1},xi::Array{Float64,1},NG::Int64,Nmid::Int64)
    # Start point
    F[1] = x[1] - xi[1]
    # first path
    for i = 2:Nmid-1
        F[i] = x[i] - C[i]*(x[i+1] - 2*x[i] + x[i-1]) - K[i]
    end
    # midpoint
    F[Nmid] = x[Nmid] - xi[2]
    # second path
    for i = Nmid+1:NG
        F[i] = x[i] - C[i]*(x[i+1] - 2*x[i] + x[i-1]) - K[i]
    end
    # end point
    F[NG+1] = x[NG+1] - xi[3]
    return(F)
end

# function to solve the system of linear equations
function linsys(x::Array{Float64,1},xprim::Array{Float64,1},λs::Array{Float64,1},ϑs::Array{Float64,1},λprim::Array{Float64,1},
                Hx::Array{SymEngine.Basic,1},Hθ::Array{SymEngine.Basic,1},Hθθ::Array{SymEngine.Basic,2},
                Hθx::Array{SymEngine.Basic,2},Δτ::Float64,NG::Int64,Nmid::Int64,H::SymEngine.Basic)
    # define relevant symbols
    X, the = symbols("X the")
    # Make array to store fixed points
    xi = fill(NaN, 3)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Array{SymEngine.Basic,1}(undef,1) # temporary hamiltonian so master isn't changed
    Hxt = Array{SymEngine.Basic,1}(undef,1)
    Hθt[1] = subs(Hθ[1], the=>0.0)
    Hxt[1] = subs(Hx[1], the=>0.0)
    Hθtt = Array{SymEngine.Basic,1}(undef,1)
    Hθttt = Array{SymEngine.Basic,1}(undef,1)
    Ht = subs(H, the=>0.0)
    Ht = subs(Ht, X=>x[Nmid]) |> float
    # loop to define fixed points
    Hxt[1] = subs(Hxt[1], X=>x[Nmid]) |> float
    Hθttt[1] = subs(Hθt[1], X=>x[Nmid]) |> float
    Hθtt[1] = subs(Hθt[1], X=>x[1]) |> float
    Hθt[1] = subs(Hθt[1], X=>x[end]) |> float
    # start point
    xi[1] = Δτ*(Hθtt[1]) + x[1]
    # midpoint
    xi[2] = -Δτ*(Hxt[1]*Ht) + x[Nmid] # x[Nmid,i] is changing between loops for some reason
    # End point
    xi[3] = Δτ*(Hθt[1]) + x[end]
    Hxθt = Array{SymEngine.Basic,2}(undef,1,1)
    # loop to calculate Hxθ
    Hxθt[1,1] = subs(Hθx[1,1], the=>0.0)
    Hxθt[1,1] = subs(Hxθt[1,1], X=>x[Nmid]) |> float
    # now transpose to get right form
    Hxθt = transpose(Hxθt) # possible source of error if ive misunderstood here
    # loop for additional midpoint terms
    xi[2] -= Δτ*(Hxθt[1,1]*Hθttt[1]) # maybe Hθttt[j]?
    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1)
    Hxt = Array{SymEngine.Basic,1}(undef,1)
    Hθθt = Array{SymEngine.Basic,2}(undef,1,1)
    Hθxt = Array{SymEngine.Basic,2}(undef,1,1)
    # Put initial values for K in
    for i = 2:NG
        K[i] = x[i]
        K[i] += Δτ*λs[i]*λprim[i]*xprim[i]
    end
    for i = 2:NG
        # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
        Hxt[1] = subs(Hx[1], X=>x[i])
        Hxt[1] = subs(Hxt[1], the=>ϑs[i]) |> float
        Hθθt[1,1] = subs(Hθθ[1,1], X=>x[i])
        Hθθt[1,1] = subs(Hθθt[1,1], the=>ϑs[i]) |> float
        Hθxt[1,1] = subs(Hθx[1,1], X=>x[i])
        Hθxt[1,1] = subs(Hθxt[1,1], the=>ϑs[i]) |> float
        # Update K's with new contributions from Hamiltonians
        K[i] -= Δτ*λs[i]*(Hθxt[1,1]*xprim[i])
        K[i] += Δτ*(Hθθt[1,1]*Hxt[1])
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
function discretise(x::Array{Float64,1},NG::Int64,Nmid::Int64)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s1 = zeros(Nmid)
    s2 = zeros(NG+2-Nmid)
    s1[1] = 0
    s2[1] = 0
    for i = 2:Nmid
        dX = x[i] - x[i-1]
        s1[i] = s1[i-1] + abs(dX)
    end
    for i = 2:(NG+2-Nmid)
        dX = x[i+Nmid-1] - x[i+Nmid-2]
        s2[i] = s2[i-1] + abs(dX)
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
    disx = zeros(NG+1)
    disx[1] = x[1]
    disx[Nmid] = x[Nmid]
    disx[NG+1] = x[NG+1]
    # This is done to linear order, which is probably good enough
    for i = 2:Nmid-1
        one = inds1[i] - 1
        two = inds1[i]
        s₀ = s1[one]
        s₁ = s1[two]
        x₀ = x[one]
        x₁ = x[two]
        disx[i] = x₀ + (ls1[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
    end

    for i = Nmid+1:NG
        one = inds2[i+1-Nmid] - 1
        two = inds2[i+1-Nmid]
        s₀ = s2[one+1-Nmid]
        s₁ = s2[two+1-Nmid]
        x₀ = x[one]
        x₁ = x[two]
        disx[i] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
    end
    return(disx)
end

# Function to generate an initial path then run the alorithm until a MAP is obtained
# This path is then returned for other functions
function gMAP(ps::Array{Float64,1},NG::Int64,Nmid::Int64,Δτ::Float64,ss1::Float64,sad::Float64,ss2::Float64)
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(ps)
    # First find the steady states and saddle point
    X1 = collect(range(ss1,stop=sad,length=Nmid))
    X2 = collect(range(sad,stop=ss2,length=Nmid))
    x = vcat(X1,X2[2:end])
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
            δ += abs(x[i] - xn[i])
        end
        println("$(δ)")
        l += 1
        if l % 200 == 0
            println("$l,$δ")
            flush(stdout)
        end
        # Now overwrite old x
        x = xn
        if δ <= 5.0*10.0^(-10)
            convrg = true
            print("$(l) steps to converge\n")
            flush(stdout)
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
    infile = "../Results/Fig3DataS/$(ARGS[1])paraS.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 4
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
    infile = "../Results/Fig3DataS/$(ARGS[1])steadS.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    l = countlines(infile)
    w = 3
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
    Δτ = 0.01
    path1 = zeros(NG+1,2)
    path2 = zeros(NG+1,2)
    # run the stability analysis for each of the hundred steady states
    for i = 1:l
        println("Run number: $(i)")
        flush(stdout)
        path1 = gMAP(ps[i,:],NG,Nmid,Δτ,steads[i,1],steads[i,2],steads[i,3])
        path2 = gMAP(ps[i,:],NG,Nmid,Δτ,steads[i,3],steads[i,2],steads[i,1])
        # define sensible names for the output files
        outfile1 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2l.csv"
        outfile2 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2h.csv"
        # open files for writing
        out_file1 = open(outfile1, "w")
        for j = 1:size(path1,1)
            line = "$(path1[j])\n"
            write(out_file1, line)
        end
        close(out_file1)
        out_file2 = open(outfile2, "w")
        for j = 1:size(path2,1)
            line = "$(path2[j])\n"
            write(out_file2, line)
        end
        close(out_file2)
    end
    return(nothing)
end

@time main()
