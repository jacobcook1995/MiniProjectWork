#!/usr/bin/env julia
# fluxaff.jl
# A script to read in the parameters and steady states and return fluxes and affinities
#
# Author: Jacob Cook
# Date: February 2019

function main()
    println("Compiled, Starting script.")
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
    
    return(nothing)
end

@time main()
