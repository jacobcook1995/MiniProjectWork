#!/usr/bin/env julia
# fsandDs.jl
# A script to read in a set of parameters for 2 species gene switch
# and the appropriate steady states and then to compare with (f1^2)/D + (f2^2)/D terms
# also should caluate the value of D at eah state
#
# Author: Jacob Cook
# Date: November 2018

# Diffusion matrix D
# ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
function D!(D::Array{Float64,2},x::Array{Float64,1},ps::Array{Float64,1})
    D[1,1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) + (ps[5]+ps[2])*x[1] + ps[6]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) + (ps[7]+ps[4])*x[2] + ps[8]
    return(D)
end

# vector of in forces
function f1!(f::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) - ps[2]*x[1]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) - ps[4]*x[2]
    return(f)
end

# vector of out forces
function f2!(f::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    f[1] = ps[5]*x[1] - ps[6]
    f[2] = ps[7]*x[2] - ps[8]
    return(f)
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
    # Only doing sum for each of three states, could expand later
    prods = zeros(l,3)
    flows = zeros(l,3)
    Ds = zeros(l,3)
    D = zeros(2,2)
    f1 = zeros(2)
    f2 = zeros(2)
    for i = 1:l
        for j = 1:3
           D = D!(D,steads[i,(2*j-1):(2*j)],ps[i,:])
           f1 = f1!(f1,steads[i,(2*j-1):(2*j)],ps[i,:])
           f2 = f2!(f2,steads[i,(2*j-1):(2*j)],ps[i,:])
           # save magnitude of each matrix
           Ds[i,j] = D[1,1]*D[2,2]
           prods[i,j] = 2*(f1[1]*f1[1]/D[1,1] + f1[2]*f1[2]/D[2,2] + f2[1]*f2[1]/D[1,1] + f2[2]*f2[2]/D[2,2])
           flows[i,j] = 4*(f1[1]*f2[1]/D[1,1] + f1[2]*f2[2]/D[2,2])
        end
    end
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])D.csv"
    out_file = open(output_file, "w")
    for i = 1:size(Ds,1)
        line = ""
        for j = 1:size(Ds,2)
            line *= "$(Ds[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])prod.csv"
    out_file = open(output_file, "w")
    for i = 1:size(prods,1)
        line = ""
        for j = 1:size(prods,2)
            line *= "$(prods[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])flow.csv"
    out_file = open(output_file, "w")
    for i = 1:size(flows,1)
        line = ""
        for j = 1:size(flows,2)
            line *= "$(flows[i,j]),"
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
