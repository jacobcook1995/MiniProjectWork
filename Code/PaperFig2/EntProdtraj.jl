#!/usr/bin/env julia
# EntProdtraj.jl
# A script to read in parameters and trajcetories and then to generate nice plots
# of the entropy productions by two differenet methods.
# 1) using the langevin entropy production term
# 2) by using a large volume reduced form of the master equation
#
# Author: Jacob Cook
# Date: January 2019

using Plots
using LaTeXStrings
using PyCall
#import PyPlot

function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig2Data/$(ARGS[1])para.csv"
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
    # now read in the 6 trajcetories
    N = 600
    traj = zeros(N+1,12)
    for i = 1:l
        for j = 1:2
            if j == 1
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])A2B.csv"
            else
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])B2A.csv"
            end
            w = 2
            open(infile, "r") do in_file
                # Use a for loop to process the rows in the input file one-by-one
                k = 1
                for line in eachline(in_file)
                    # parse line by finding commas
                    L = length(line)
                    comma = fill(0,w+1)
                    m = 1
                    for l = 1:L
                        if line[l] == ','
                            m += 1
                            comma[m] = l
                        end
                    end
                    comma[end] = L+1
                    for l = 1:w
                        traj[k,4*(i-1)+2*(j-1)+l] = parse(Float64,line[(comma[l]+1):(comma[l+1]-1)])
                    end
                    k += 1
                end
            end
        end
    end
    
    return(nothing)
end

@time main()
