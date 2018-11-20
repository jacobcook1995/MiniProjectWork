#!/usr/bin/env julia
# entgill.jl
# A script to read in parameters and steady states and then do Gillespie simulationa
# in order to obtain values of entropy for each steady state
#
# Author: Jacob Cook
# Date: November 2018

using SparseArrays

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars[1] += 1 # A produced
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
    elseif r < sum(rs[1:3])
        vars[1] -= 1 # A decays
    elseif r < sum(rs[1:4])
        vars[1] += 1 # A regenerated
    elseif r < sum(rs[1:5])
        vars[2] += 1 # B produced
    elseif r < sum(rs[1:6])
        vars[2] -= 1 # B unravels
    elseif r < sum(rs[1:7])
        vars[2] -= 1 # B decays
    else
        vars[2] += 1 # B regenerated
    end
    return(vars)
end

function gillespie(stead::Array{Float64,1},mid::Array{Float64,1},ps::Array{Float64,1},noits::Int64,Ω::Int64)
    # extract all parameters and then rescale
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    k = ps[1]*Ω
    kmin = ps[2]
    q = ps[3]*Ω
    qmin = ps[4]
    K = ps[5]
    Kmin = ps[6]*Ω
    Q = ps[7]
    Qmin = ps[8]*Ω
    r = ps[9]
    f = ps[10]/(Ω^2)
    # make relevant data structures
    times = zeros(2)
    vars = fill(0,2,2)
    sad = round.(Int64,mid*Ω)
    star = round.(Int64,stead*Ω)
    vars[:,1] = star
    hist = spzeros(80*Ω,80*Ω) # Possible source of error
    # run gillespie loop
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        hist[vars[1,1]+1,vars[2,1]+1] += times[2] - times[1]
        # do gillepsie step
        vars[:,2] = step(rs,vars[:,1])
        vars[:,1] = vars[:,2]
        times[1] = times[2]
    end
    # normalise probability distribution
    hist = hist/times[2]
    # Now extract entropy of the state from the histogram
    SA = 0.0
    SB = 0.0
    # rescale for first state
    histn = hist[1:(sad[1]+1),(sad[2]+1):end]/sum(hist[1:(sad[1]+1),(sad[2]+1):end])
    for j = 1:size(histn,2)
        for i = 1:size(histn,1)
            if histn[i,j] != 0 && ~isnan(histn[i,j])
                SB -= histn[i,j]*log(histn[i,j])
            end
        end
    end
    # rescale for second state
    histn = hist[(sad[1]+1):end,1:(sad[2]+1)]/sum(hist[(sad[1]+1):end,1:(sad[2]+1)])
    for j = 1:size(histn,2)
        for i = 1:size(histn,1)
            if histn[i,j] != 0 && ~isnan(histn[i,j])
                SA -= histn[i,j]*log(histn[i,j])
            end
        end
    end
    # return values in appropriate order
    if stead[1] > mid[1]
        return(SA,SB)
    else
        return(SB,SA)
    end
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
    S = zeros(l,2)
    # parameters to fiddle for the gillespie
    Ω = 1
    noits = 1000000000
    for i = 1:l
        println(i)
        flush(stdout)
        S[i,1], S[i,2] = gillespie([steads[i,1],steads[i,2]],[steads[i,3],steads[i,4]],ps[i,:],noits,Ω)
    end
    # Now write out ents to file
    output_file = "../Results/Fig3Data/$(ARGS[1])ent.csv"
    out_file = open(output_file, "w")
    for i = 1:size(S,1)
        line = ""
        for j = 1:size(S,2)
            line *= "$(S[i,j]),"
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
