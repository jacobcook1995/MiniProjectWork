#!/usr/bin/env julia
# fluxaff.jl
# A script to read in the parameters and steady states and return fluxes and affinities
#
# Author: Jacob Cook
# Date: February 2019

# function to construct the rates
function rates(A::Float64,B::Float64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*B), kmin*A, K*A, Kmin, r*q/(r + f*A*A), qmin*B, Q*B, Qmin ]
    return(rates)
end

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
    # first find all fluxes and affinities of the states
    fluxes = zeros(l,8)
    affs = zeros(l,8)
    for i = 1:l
        k = ps[i,1]
        kmin = ps[i,2]
        q = ps[i,3]
        qmin = ps[i,4]
        K = ps[i,5]
        Kmin = ps[i,6]
        Q = ps[i,7]
        Qmin = ps[i,8]
        r = ps[i,9]
        f = ps[i,10]
        # First steady state
        rs = rates(steads[i,1],steads[i,2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        fluxes[i,1] = rs[1] - rs[2]
        fluxes[i,2] = rs[3] - rs[4]
        fluxes[i,3] = rs[5] - rs[6]
        fluxes[i,4] = rs[7] - rs[8]
        affs[i,1] = log(rs[1]/rs[2])
        affs[i,2] = log(rs[3]/rs[4])
        affs[i,3] = log(rs[5]/rs[6])
        affs[i,4] = log(rs[7]/rs[8])
        # Second steady state
        rs = rates(steads[i,5],steads[i,6],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        fluxes[i,5] = rs[1] - rs[2]
        fluxes[i,6] = rs[3] - rs[4]
        fluxes[i,7] = rs[5] - rs[6]
        fluxes[i,8] = rs[7] - rs[8]
        affs[i,5] = log(rs[1]/rs[2])
        affs[i,6] = log(rs[3]/rs[4])
        affs[i,7] = log(rs[5]/rs[6])
        affs[i,8] = log(rs[7]/rs[8])
    end
    # Literally just going to output everything
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])fluxes.csv"
    out_file = open(output_file, "w")
    for i = 1:size(fluxes,1)
        line = ""
        for j = 1:size(fluxes,2)
            line *= "$(fluxes[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    close(out_file)
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])affs.csv"
    out_file = open(output_file, "w")
    for i = 1:size(affs,1)
        line = ""
        for j = 1:size(affs,2)
            line *= "$(affs[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    close(out_file)
    return(nothing)
end

@time main()
