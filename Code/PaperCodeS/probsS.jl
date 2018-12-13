#!/usr/bin/env julia
# probsS.jl
# A script to calculate the probability of being in either steday state via gillespie simulation
# this is then output and saved as a .csv file
# This is now done for the 1 species Schlögl model
#
# Author: Jacob Cook
# Date: December 2018

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
    elseif r < rs[1] + rs[2] + rs[3]
        vars -= 1 # X decays
    else
        vars += 1 # X regenerated
    end
    return(vars)
end

function gillespie(stead::Float64,mid::Float64,ps::Array{Float64,1},noits::Int64,Ω::Float64)
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
    th = 0
    tl = 0
    if star > sad
        high = true
    else
        high = false
    end
    N = 0 # number of switches
    # run gillespie loop
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1],k1,K1,k2,K2)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        if high == true
            th += times[2] - times[1]
        else
            tl += times[2] - times[1]
        end
        # do gillepsie step
        vars[2] = step(rs,vars[1])
        # now need step to determine if state has switched
        if high == true && vars[2] <= sad
            high = false
            N += 1
        elseif high == false && vars[2] >= sad
            high = true
            N += 1
        end
        vars[1] = vars[2]
        times[1] = times[2]
    end
    # find probabilities from times
    T = tl+th
    ph = th/T
    pl = tl/T
    return(ph,pl,N,T)
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
    # now run gillespie to find plow and phigh
    probs = zeros(l,2)
    Ωs = zeros(l)
    noits = 100000000
    thresh = 5000
    for i = 1:l
        # need a step here to ensure that a resonable volume system is simulated
        println("Run $(i)")
        flush(stdout)
        tot = (steads[i,1] + steads[i,3])/2
        println(tot)
        # this step is somwhat provisional, NEED TO THINK MORE
        if tot >= 50.0
            Ωs[i] = 0.125
        elseif tot >= 25.0
            Ωs[i] = 0.25
        elseif tot >= 10.0
            Ωs[i] = 0.5
        elseif tot >= 7.5
            Ωs[i] = 0.75
        elseif tot >= 5.0
            Ωs[i] = 1.0
        elseif tot >= 2.0
            Ωs[i] = 2.0
        else
            Ωs[i] = 3.0
        end
        repeat = false
        probs[i,1], probs[i,2], N, T = gillespie(steads[i,1],steads[i,2],ps[i,:],noits,Ωs[i])
        if N < thresh
            repeat = true
        end
        count = 0
        while repeat == true
            if count % 2 == 0
                pht, plt, n, τ = gillespie(steads[i,3],steads[i,2],ps[i,:],noits,Ωs[i])
            else
                pht, plt, n, τ = gillespie(steads[i,1],steads[i,2],ps[i,:],noits,Ωs[i])
            end
            count += 1
            probs[i,1] = (probs[i,1]*T + pht*τ)/(T+τ)
            probs[i,2] = (probs[i,2]*T + plt*τ)/(T+τ)
            T += τ
            N += n
            if N > thresh
                repeat = false
            end
        end
    end
    # now output data
    outfile = "../Results/Fig3DataS/$(ARGS[1])probsS.csv"
    out_file = open(outfile, "w")
    for i = 1:l
        line = "$(probs[i,1]),$(probs[i,2]),$(Ωs[i])\n"
        write(out_file,line)
    end
    return(nothing)
end

@time main()
