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

function gillespie(stead1::Float64,mid::Float64,stead2::Float64,ps::Array{Float64,1},noits::Int64,Ω::Float64)
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
    star = round(Int64,stead1*Ω)
    fin = round(Int64,stead2*Ω)
    vars[1] = star
    th = 0
    tl = 0
    if star >= sad
        high = true
        if sad - fin < 2
            sad = sad + 2
        end
    else
        high = false
        if sad - star < 2
            sad = sad + 2 # step to avoid overlapping saddle points
        end
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
        if high == true && vars[2] <= sad # need to establish what is reasonable here
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
    # read in actions
    acts = zeros(l,8)
    datayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2lD.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2hD.csv"
            end
            # check it exits
            if isfile(infile)
                w = 4
                open(infile, "r") do in_file
                    # only one line now
                    for line in eachline(in_file)
                        # parse line by finding commas
                        L = length(line)
                        comma = fill(0,w+1)
                        m = 1
                        for k = 1:L
                            if line[k] == ','
                                m += 1
                                comma[m] = k
                            end
                        end
                        comma[end] = L+1
                        for k = 1:w
                            acts[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
                        end
                    end
                end
            else # if file doesn't exist inform user of missing data and exclude from plots
                if j == 1
                    println("Missing data for run $(i) high A to high B")
                    datayn[i] = false
                else
                    println("Missing data for run $(i) high B to high A")
                    datayn[i] = false
                end
            end
        end
    end
    # now run gillespie to find plow and phigh
    probs = zeros(l,2)
    Ωs = zeros(l)
    noits = 100000000
    thresh = 50000
    # These comments are obviously provisonal as they might run into slower or faster patches still
    inds = 1:17 # Reached 7 # done 6 out of 17
    inds = 18:34 # Reached 23 # done 5 out of 17
    inds = 35:51 # Reached 41 # done 6 out of 17
    inds = 52:68 # Reached 59 # done 7 out of 17
    inds = 69:85 # Reached 73 # done 4 out of 17
    inds = 86:100 # Reached 92 # done 6 out of 15
    for i = inds
        # need a step here to ensure that a resonable volume system is simulated
        println("Run $(i)")
        flush(stdout)
        maxA = maximum([acts[i,2],acts[i,6]])
        if maxA < 0.05
            Ωs[i] = 7.0/maxA
        elseif maxA > 2.5
            Ωs[i] = 11.0/maxA
        else
            Ωs[i] = 10.0/maxA
        end
        # maybe looking at far too small volumes here
        repeat = false
        probs[i,1], probs[i,2], N, T = gillespie(steads[i,1],steads[i,2],steads[i,3],ps[i,:],noits,Ωs[i])
        println("N = $(N)")
        flush(stdout)
        if N < thresh
            repeat = true
        end
        count = 0
        while repeat == true
            if count % 2 == 0
                pht, plt, n, τ = gillespie(steads[i,3],steads[i,2],steads[i,1],ps[i,:],noits,Ωs[i])
            else
                pht, plt, n, τ = gillespie(steads[i,1],steads[i,2],steads[i,3],ps[i,:],noits,Ωs[i])
            end
            count += 1
            probs[i,1] = (probs[i,1]*T + pht*τ)/(T+τ)
            probs[i,2] = (probs[i,2]*T + plt*τ)/(T+τ)
            T += τ
            N += n
            println("N = $(N)")
            flush(stdout)
            if N > thresh
                repeat = false
            end
        end
    end
    # now output data
    outfile = "../Results/Fig3DataS/$(ARGS[1])probsS.csv"
    # Should have a step here so that I can replace lines
    temprobs = zeros(l,3)
    if isfile(outfile)
        # now read in steady states
        l = countlines(outfile)
        w = 3
        steads = zeros(l,w)
        open(outfile, "r") do in_file
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
                    temprobs[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
    end
    for i = 1:l
        if ~any(x->x==i, inds)
            probs[i,:] = temprobs[i,1:2]
            Ωs[i] = temprobs[i,3]
        end
    end
    out_file = open(outfile, "w")
    for i = 1:l
        line = "$(probs[i,1]),$(probs[i,2]),$(Ωs[i])\n"
        write(out_file,line)
    end
    return(nothing)
end

@time main()
