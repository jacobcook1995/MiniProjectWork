#!/usr/bin/env julia
# steadstates.jl
# A script to read in a set of parameters for 2 species gene switch
# then use these parameters to find steady states and saddle points
# Finally these steady states should be output as a csv file
#
# Author: Jacob Cook
# Date: November 2018

using Roots

# function to test if parameter set allows three +ve stationary points
function nullcline(ps::Array{Float64,1})
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    # A1(x) = sqrt((r/f)*(q/((qmin+Q)*B - Qmin) - 1))
    A1(x) = real(sqrt(complex((ps[9]/ps[10])*(ps[3]/((ps[4]+ps[7])*x - ps[8]) - 1))))
    # A2(x) = (1/(kmin+K))*(kr/(r+f*x^2) + Kmin)
    A2(x) = (1/(ps[2]+ps[5]))*((ps[1]*ps[9])/(ps[9]+ps[10]*x^2) + ps[6])
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    q = ps[3]
    Q = ps[7]
    while three == false
        bs1 = fzeros(g, 0.0, 0.1) # catches zeros about the origin
        bs2 = fzeros(g, 0.1, 2*q/Q) # the upper bound here is potentially problematic
        bs = vcat(bs1,bs2)
        # bs = fzeros(g, 0, 15.0)
        n = length(bs)
        gs = 0
        bad = zeros(Int64,0)
        for i = 1:n
            # check if this is actual solution and not artifact of the numerical method
            gs = g(bs[i])
            tol = 1.0e-14
            if gs >= 0 + tol || gs <= 0 - tol
                bad = append!(bad,i)
            end
        end
        if length(bad) != 0
            n = n - length(bad)
            bs = deleteat!(bs, bad)
        end
        if n == 3
            three = true
        end
    end
    # vector of steady states
    sss = [ A1(bs[1]), bs[1], A1(bs[2]), bs[2], A1(bs[3]), bs[3] ]
    # define three steady states as a single vector
    return(sss)
end

# function to find schnakenberg entropy productions
# ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
function sch(ps::Array{Float64,1},A::Float64,B::Float64)
    r1 = ps[1]*ps[9]/(ps[9]+ps[10]*B^2)
    r2 = ps[5]*A
    r3 = ps[3]*ps[9]/(ps[9]+ps[10]*A^2)
    r4 = ps[7]*B
    rmin1 = ps[2]*A
    rmin2 = ps[6]
    rmin3 = ps[4]*B
    rmin4 = ps[8]
    F1 = r1 - rmin1
    F2 = r2 - rmin2
    F3 = r3 - rmin3
    F4 = r4 - rmin4
    A1 = log(r1/rmin1)
    A2 = log(r2/rmin2)
    A3 = log(r3/rmin3)
    A4 = log(r4/rmin4)
    S1 = F1*A1
    S2 = F2*A2
    S3 = F3*A3
    S4 = F4*A4
    S12 = S1 + S2
    S34 = S3 + S4
    S = S12 + S34
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
    # find steady states
    steads = zeros(l,6)
    for i = 1:l
        steads[i,:] = nullcline(ps[i,:])
    end
    # find entropies
    entps = zeros(l,3)
    for i = 1:l
        entps[i,1] = sch(ps[i,:],steads[i,1],steads[i,2])
        entps[i,2] = sch(ps[i,:],steads[i,3],steads[i,4])
        entps[i,3] = sch(ps[i,:],steads[i,5],steads[i,6])
    end
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])stead.csv"
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
    output_file = "../Results/Fig3Data/$(ARGS[1])schnak.csv"
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
