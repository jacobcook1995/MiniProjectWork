#!/usr/bin/env julia
# steadstatesS.jl
# A script to read in a set of parameters for 1 species Schlögl model
# then use these parameters to find steady states and saddle points
# Finally these steady states should be output as a .csv file
#
# Author: Jacob Cook
# Date: December 2018

using Roots

# function to test if parameter set allows three +ve stationary points
function nullcline(ps::Array{Float64,1})
    # ps = [ k1, K1, k2, K2 ]
    a = ps[3]
    b = -ps[4]
    c = ps[2]
    d = -ps[1]
    f(x) = a*x^3 + b*x^2 + c*x + d
    Xs = fzeros(f, 0.0, ps[4]/ps[3] + ps[1]/ps[2])
    # vector of steady states
    sss = [ Xs[3], Xs[2], Xs[1] ]
    # define three steady states as a single vector
    return(sss)
end

# function to find schnakenberg entropy productions
# ps = [ k1, K1, k2, K2 ]
function sch(ps::Array{Float64,1},X::Float64)
    r1 = ps[1]
    r2 = ps[3]*X^3
    rmin1 = ps[2]*X
    rmin2 = ps[4]*X^2
    F1 = r1 - rmin1
    F2 = r2 - rmin2
    A1 = log(r1/rmin1)
    A2 = log(r2/rmin2)
    S1 = F1*A1
    S2 = F2*A2
    S = S1 + S2
    return(S)
end

function main()
    println("Compiled, Starting script.")
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
    # find steady states
    steads = zeros(l,3)
    for i = 1:l
        steads[i,:] = nullcline(ps[i,:])
    end
    # find entropies
    entps = zeros(l,3)
    for i = 1:l
        entps[i,1] = sch(ps[i,:],steads[i,1])
        entps[i,2] = sch(ps[i,:],steads[i,2])
        entps[i,3] = sch(ps[i,:],steads[i,3])
    end
    # Finally write out to file
    output_file = "../Results/Fig3DataS/$(ARGS[1])steadS.csv"
    out_file = open(output_file, "w")
    for i = 1:size(steads,1)
        line = ""
        for j = 1:size(steads,2)
            line *= "$(steads[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    # Finally write out to file
    output_file = "../Results/Fig3DataS/$(ARGS[1])schnakS.csv"
    out_file = open(output_file, "w")
    for i = 1:size(entps,1)
        line = ""
        for j = 1:size(entps,2)
            line *= "$(entps[i,j]),"
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
