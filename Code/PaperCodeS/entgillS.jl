#!/usr/bin/env julia
# entgillS.jl
# A script to read in parameters and steady states and then do Gillespie simulationa
# in order to obtain values of entropy production based on master equation dynamics
#
# Author: Jacob Cook
# Date: December 2018

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64)
    rates = [ k1, K1*X, k2*X*(X-1)*(X-2), K2*X*(X-1) ]
    return(rates)
end
# function to construct the rates
function rates(X::Float64,k1::Float64,K1::Float64,k2::Float64,K2::Float64)
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
function step(rates::Array{Float64,1},vars::Int64,reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars += 1 # X produced
        p = rs[1]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels
        p = rs[2]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars -= 1 # X decays
        p = rs[3]
        reac = 3
    else
        vars += 1 # X regenerated
        p = rs[4]
        reac = 4
    end
    return(vars,p,reac)
end

# function to find reverse probability
function rev(rs::Array{Float64,1},reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 5
        if reac % 2 == 0 # even rates backwards
            return(rs[reac-1])
        else # odd forward
            return(rs[reac+1])
        end
    else
        error("Invalid reaction code returned")
    end
end

function gillespie(steadXh::Float64,steadXl::Float64,mid::Float64,ps::Array{Float64,1},noits::Int64,Ω::Int64,pf::Array{Float64,1},pb::Array{Float64,1})
    # extract all parameters and then rescale
    # ps = [ k1, K1, k2, K2 ]
    k1 = ps[1]*Ω
    K1 = ps[2]
    k2 = ps[3]/(Ω^2)
    K2 = ps[4]/Ω
    # set up these gloablly
    times = zeros(2)
    vars = fill(0,2)
    Sh = Sl = 0
    # setup for finding if left steady state
    sad = round(Int64,mid*Ω)
    fin = false
    while fin == false
        crossX = false
        times[1] = 0
        vars[1] = round(Int64,steadXh*Ω)
        reac = 0
        minX = vars[1]
        maxX = vars[1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1],k1,K1,k2,K2)
            if i != 1
                # now use reac to calculate reverse rate
                pb[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[2], pf[i], reac = step(rs,vars[1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[2],k1,K1,k2,K2)
                pb[end] = rev(rs,reac)
            end
            # check if max and minimums should be updated
            if vars[2] < minX
                minX = vars[2]
            elseif vars[2] > maxX
                maxX = vars[2]
            end
            vars[1] = vars[2]
            times[1] = times[2]
        end
        Sh = 0
        for i = 1:length(pb)
            Sh += log(pf[i]) - log(pb[i])
        end
        Sh = Sh/(Ω*times[2])
        # check simulated trajectories remain in bounds
        if (steadXh - mid) < 0
            if sad < maxX
                crossX = true
            end
        else
            if sad > minX
                crossX = true
            end
        end
        if crossX == false
            fin = true
            break
        else
            println(steadXh)
            println(mid)
            flush(stdout)
        end
    end

    # Second run for other state
    fin = false
    while fin == false
        crossX = false
        times[1] = 0
        vars[1] = round(Int64,steadXl*Ω)
        reac = 0
        minX = vars[1]
        maxX = vars[1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1],k1,K1,k2,K2)
            if i != 1
                # now use reac to calculate reverse rate
                pb[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[2], pf[i], reac = step(rs,vars[1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[2],k1,K1,k2,K2)
                pb[end] = rev(rs,reac)
            end
            # check if max and minimums should be updated
            if vars[2] < minX
                minX = vars[2]
            elseif vars[2] > maxX
                maxX = vars[2]
            end
            vars[1] = vars[2]
            times[1] = times[2]
        end
        Sl = 0
        for i = 1:length(pb)
            Sl += log(pf[i]) - log(pb[i])
        end
        Sl = Sl/(Ω*times[2])
        # check simulated trajectories remain in bounds
        if (steadXl - mid) < 0
            if sad < maxX
                crossX = true
            end
        else
            if sad > minX
                crossX = true
            end
        end
        if crossX == false
            fin = true
            break
        else
            println(steadXl)
            println(mid)
            flush(stdout)
        end
    end
    return([Sh,Sl])
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
    ents = zeros(l,2)
    # things to fiddle with to get to work
    noits = 100000000
    Ω = 5000
    pf = zeros(noits)
    pb = zeros(noits)
    for i = 1:l
        println(i)
        flush(stdout)
        ents[i,:] = gillespie(steads[i,1],steads[i,3],steads[i,2],ps[i,:],noits,Ω,pf,pb)
    end
    # Now write out ents to file
    output_file = "../Results/Fig3DataS/$(ARGS[1])gillS.csv"
    out_file = open(output_file, "w")
    for i = 1:size(ents,1)
        line = ""
        for j = 1:size(ents,2)
            line *= "$(ents[i,j]),"
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
