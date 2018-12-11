#!/usr/bin/env julia
# steadentS.jl
# A script to read in parameters and steady states and then do Gillespie simulationa
# in order to obtain values of entropy for each steady state
#
# Author: Jacob Cook
# Date: December 2018

using SparseArrays

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64)
    rates = [ k1, K1*X, k2*X*(X-1)*(X-2), K2*X*(X-1) ]
    return(rates)
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Int64)
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars += 1 # X produced
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels
    elseif r < sum(rs[1:3])
        vars -= 1 # X decays
    else
        vars += 1 # X regenerated
    end
    return(vars)
end

function gillespie(stead::Float64,mid::Float64,ps::Array{Float64,1},noits::Int64,Ω::Int64)
    # extract all parameters and then rescale
    # ps = [ k1, K1, k2, K2 ]
    k1 = ps[1]*Ω
    K1 = ps[2]
    k2 = ps[3]/(Ω^2)
    K2 = ps[4]/Ω
    # make relevant data structures
    times = zeros(2)
    vars = fill(0,2)
    sad = round(Int64,mid*Ω)
    star = round(Int64,stead*Ω)
    vars[1] = star
    hist = spzeros(250*Ω) # making matrix huge as now using a sparse array
    # run gillespie loop
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1],k1,K1,k2,K2)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        hist[vars[1]+1] += times[2] - times[1]
        # do gillepsie step
        vars[2] = step(rs,vars[1])
        vars[1] = vars[2]
        times[1] = times[2]
    end
    # normalise probability distribution
    hist = hist/times[2]
    # Now extract entropy of the state from the histogram
    Sh = 0.0
    Sl = 0.0
    # rescale for first state
    histn = hist[1:(sad+1)]/sum(hist[1:(sad+1)])
    for i = 1:size(histn,1)
        if histn[i] != 0 && ~isnan(histn[i])
            Sl -= histn[i]*log(histn[i])
        end
    end
    # rescale for second state
    histn = hist[(sad+1):end]/sum(hist[(sad+1):end])
    for i = 1:size(histn,1)
        if histn[i] != 0 && ~isnan(histn[i])
            Sh -= histn[i]*log(histn[i])
        end
    end
    return(Sh,Sl)
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
    S = zeros(l,2)
    # parameters to fiddle for the gillespie
    Ω = 1
    noits = 1000000000
    for i = 1:l
        println(i)
        flush(stdout)
        S[i,1], S[i,2] = gillespie(steads[i,1],steads[i,2],ps[i,:],noits,Ω)
        # step to catch 0s
        if S[i,1] == 0.0
            println("S1 rerun")
            _, S[i,1] = gillespie(steads[i,3],steads[i,2],ps[i,:],noits,Ω)
        elseif S[i,2] == 0.0
            println("S2 rerun")
            S[i,2], _ = gillespie(steads[i,3],steads[i,2],ps[i,:],noits,Ω)
        end
    end
    # Now write out ents to file
    output_file = "../Results/Fig3DataS/$(ARGS[1])entS.csv"
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
