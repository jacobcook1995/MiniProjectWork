#!/usr/bin/env julia
# gill1spec.jl
# A script that performs a gillespie simulation of the 1 species schlögl model
# with explicit promotor un/binding. Histograms of the occupation and the waiting time
# distributions of the two states are generated and output
#
# Author: Jacob Cook
# Date: December 2018

using Roots
using Plots
import GR # this is necessary to avoid a world age error when using GR in function

# function to provide normalised parameters for the gillespie simulation
function paras(Ω::Int64)
    # just gonna hard code some in for the moment
    k1 = 4.26
    K1 = 9.85
    k2 = 0.115
    K2 = 3.04
    # then normalise appropriatly
    k1 = k1*Ω
    K1 = K1
    k2 = k2/(Ω^2)
    K2 = K2/Ω
    return(k1,K1,k2,K2)
end

# function to generate steady states
function states(k1::Float64,K1::Float64,k2::Float64,K2::Float64,Ω::Int64)
    # unnormalise the parameters as this is at Langevin level
    k1 = k1/Ω
    K1 = K1
    k2 = k2*(Ω^2)
    K2 = K2*Ω
    a = k2
    b = -K2
    c = K1
    d = -k1
    f(x) = a*x^3 + b*x^2 + c*x + d
    Xs = fzeros(f, 0.0, k1/K1 + K2/k2)
    # vector of steady states
    sss = [ Xs[3], Xs[2], Xs[1] ]
    # define three steady states as a single vector
    println(sss)
    flush(stdout) # needed to ensure output prints when using nohup
    return(sss)
end


# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,diffs::Bool=false)
    rates = [ k1, K1*X, k2*X*(X-1)*(X-2), K2*X*(X-1) ]
    if diffs == false
        return(rates)
    else
        dX = [ 1, -1, -1, 1 ]
        return(rates,dX)
    end
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step!(rates::Array{Float64,1},vars::Int64)
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
    return(vars) # there's always the potential of inverting this to speed the code up
end

# function to actually run a gillepsie simulation
function gillespie(k1::Float64,K1::Float64,k2::Float64,K2::Float64,noits::Int64,star::Int64,Ω::Int64,fin::Int64)
    # define high X and lowX states for comparison
    hX = star
    lX = fin
    # set bool as being in high state initially
    high = true
    # make vector to store waiting times
    wh = Array{Float64,1}(undef,0)
    wl = Array{Float64,1}(undef,0)
    tp = 0 # prior time
    times = zeros(2)
    vars = fill(0,2)
    vars[2] = star
    hist = zeros(100*Ω)
    for i = 1:noits
        vars[1] = vars[2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[1],k1,K1,k2,K2)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[2] = step!(rs,vars[1])
        # add to histogram
        hist[vars[1]+1] += times[2] - times[1]
        # now step to check if transistion has occurred
        if vars[2] >= hX && high == false
            high = true
            wl = vcat(wl,times[2]-tp)
            tp = times[2]
        elseif vars[2] <= lX && high == true
            high = false
            wh = vcat(wh,times[2]-tp)
            tp = times[2]
        end
    end
    hist = hist/times[2]
    return(hist,wh,wl,times[2])
end

# function to printout my data
function printout(hist::Array{Float64,1},wh::Array{Float64,1},wl::Array{Float64,1},Time::Float64)
    if isfile("../Results/wh$(ARGS[1])V$(ARGS[2])S.csv")
        # if file already exists should add to bottom of it
        output_file = "../Results/wh$(ARGS[1])V$(ARGS[2])S.csv"
        out_file = open(output_file, "a")
        for i = 1:size(wh,1)
            line = "$(wh[i])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        output_file = "../Results/wh$(ARGS[1])V$(ARGS[2])S.csv"
        out_file = open(output_file, "w")
        for i = 1:size(wh,1)
            line = "$(wh[i])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    if isfile("../Results/wl$(ARGS[1])V$(ARGS[2])S.csv")
        # if file already exists should add to bottom of it
        output_file = "../Results/wl$(ARGS[1])V$(ARGS[2])S.csv"
        out_file = open(output_file, "a")
        for i = 1:size(wl,1)
            line = "$(wl[i])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        output_file = "../Results/wl$(ARGS[1])V$(ARGS[2])S.csv"
        out_file = open(output_file, "w")
        for i = 1:size(wl,1)
            line = "$(wl[i])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    # histogram in this case is one dimensional so not hard to represent
    if isfile("../Results/hist$(ARGS[1])V$(ARGS[2])S.csv")
        # need to read in old file and then add new histogram to it and then output
        input_file = "../Results/hist$(ARGS[1])V$(ARGS[2])S.csv"
        vol = convert(Int64,(countlines(input_file)) - 3)
        histr = Array{Float64,1}(undef,vol)
        a = 0
        i = 1
        T = 0
        open(input_file, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                if line[1] == '#'
                    a += 1
                    i = 1
                else
                    if a == 2
                        T = parse(Float64,line)
                    else
                        histr[i] = parse(Float64,line)
                        i += 1
                    end
                end
            end
        end
        hist = (Time*hist .+ T*histr)/(T + Time)
        # update time to include time of new run
        Time += T
    end
    # if file doesn't exist should make it
    output_file = "../Results/hist$(ARGS[1])V$(ARGS[2])S.csv"
    out_file = open(output_file, "w")
    line = "# Start #\n"
    write(out_file,line)
    for i = 1:size(hist,1)
        line = "$(hist[i])\n"
        write(out_file,line)
    end
    line = "# Time #\n"
    write(out_file,line)
    line = "$(Time)\n"
    write(out_file,line)
    close(out_file)
    return(nothing)
end

function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    # then check that a system volume has been provided
    elseif length(ARGS) == 1
        println("Error: Need to provide an argument to set system volume.")
        return(nothing)
    end
    # Then take system volume Ω, check if provided value is integer
    Ω = 0
    try Ω = parse(Int64,ARGS[2])
    catch y
        if isa(y, ArgumentError) # would only really expect an argument error
            println("Error: System volume has to be integer.")
            return(nothing)
        end
    end
    # now want to obtain the parameters of the simulation
    k1, K1, k2, K2 = paras(Ω)
    # use these parameters to generate the steady states
    sss = states(k1,K1,k2,K2,Ω)
    # scale starting posisition by volume
    star = round(Int64,sss[1]*Ω)
    fin = round(Int64,sss[3]*Ω)
    # Now ready to set the gillespie simulation running
    noits = 100000000#00
    hist, wh, wl, T = gillespie(k1,K1,k2,K2,noits,star,Ω,fin)
    # now printout
    printout(hist,wh,wl,T)
    ps = [ k1, K1, k2, K2, Ω ]
    output_file = "../Results/ps$(ARGS[1])V$(ARGS[2])S.csv"
    out_file = open(output_file, "w")
    for i = 1:length(ps)
        line = "$(ps[i])\n"
        write(out_file,line)
    end
    close(out_file)
    return(nothing)
end

@time main()
