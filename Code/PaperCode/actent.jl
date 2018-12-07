#!/usr/bin/env julia
# actent.jl
# A script to read in the parameters and paths and then use to determine
# entropy produced along the path and action of forward and reverse paths
#
# Author: Jacob Cook
# Date: December 2018

using SymEngine

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
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
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
function gensyms(ps::Array{Float64,1})
    # create symbolic objects
    ϑ = ϑs()
    λ = λs()
    b = bs()
    Dmin = Dmins()
    # specify symbols that will be substituted for
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("K k Q q kmin Kmin qmin Qmin f r")
    # now perform substitutions
    λ = subs(λ, k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
    λ = subs(λ, Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    for i = 1:2
        ϑ[i] = subs(ϑ[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        ϑ[i] = subs(ϑ[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        b[i] = subs(b[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        b[i] = subs(b[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    end
    for j = 1:2
        for i = 1:2
            Dmin[i,j] = subs(Dmin[i,j], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
            Dmin[i,j] = subs(Dmin[i,j], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        end
    end
    return(ϑ,λ,b,Dmin)
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
    # 2 for loops to loop over every file
    for i = 1#:l
        for j = 1:2
            if j == 1
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2B.csv"
            else
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2A.csv"
            end
            # check if file
            if isfile(infile)
                # now should read in path
                l = countlines(infile)
                w = 2
                path = zeros(l,w)
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
                            path[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                        end
                        k += 1
                    end
                end

            else # Tell users about missing files
                if j == 1
                    println("No A to B path for $(i)")
                else
                    println("No B to A path for $(i)")
                end
            end
        end
    end
    return(nothing)
end

@time main()
