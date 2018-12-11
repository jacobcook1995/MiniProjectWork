#!/usr/bin/env julia
# parasetS.jl
# A script to generate and then output a set of parameters for 1 species schlögl
# these parameters are then output as a .csv file
#
# Author: Jacob Cook
# Date: December 2018

using Roots

# function to generate unifrom random numbers between from and to
function myrand(to::Float64,from::Float64)
    r = rand()*(to-from) + from
    return(r)
end

# Need to think of a good parameter randomisation scheme
function paras()
    works = false
    paras = zeros(4)
    while works == false
        # k randomly selected between 1.0 and 100.0
        k1 = myrand(10.0,0.1)
        K1 = myrand(10.0,0.1)
        # q then smaller than k by a ratio between 1.0 and 0.01
        k2 = myrand(10.0,0.1)
        K2 = myrand(10.0,0.1)
        # now gather parameters
        paras = [ k1, K1, k2, K2 ]
        works = steady3(paras)
    end
    return(paras)
end

# function to test if parameter set allows three +ve stationary points
function steady3(ps::Array{Float64,1})
    # ps = [ k1, K1, k2, K2 ]
    # find determinent of the cubic equation
    a = ps[3]
    b = -ps[4]
    c = ps[2]
    d = -ps[1]
    Δ = 18*a*b*c*d - 4*(b^3)*d + (b^2)*(c^2) - 4*a*(c^3) - 27*(a^2)*(d^2)
    if Δ > 0.0
        return(true)
    else
        return(false)
    end
end

function main()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    # then check that a system volume has been provided
    elseif length(ARGS) == 1
        println("Error: Need to provide an argument to set number of parameter sets.")
        return(nothing)
    end
    # Then take number of parameters sets N, check if provided value is integer
    N = 0
    try N = parse(Int64,ARGS[2])
    catch y
        if isa(y, ArgumentError) # would only really expect an argument error
            println("Error: Number of parameter sets has to be integer.")
            return(nothing)
        end
    end
    # ps = [ k1, K1, k2, K2 ]
    ps = zeros(N,4) # matrix to store full parameter set
    for i = 1:N
        println("Set: $(i)")
        ps[i,:] = paras()
    end
    # Finally write out to file
    output_file = "../Results/Fig3DataS/$(ARGS[1])paraS.csv"
    out_file = open(output_file, "w")
    for i = 1:size(ps,1)
        line = ""
        for j = 1:size(ps,2)
            line *= "$(ps[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    return(nothing)
end

@time main()
