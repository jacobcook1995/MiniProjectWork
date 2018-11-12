#!/usr/bin/env julia
# trajentsch.jl
# A script to take in a trajectory from the Langevin formulalism convert to a sequence
# of jumps in the master equation formulation, and then calculate the entropy production of this process
# This is now done for the Schlögl model as it's easier to interpret
#
# Author: Jacob Cook
# Date: November 2018
using Plots

# function that takes in two points and finds the probability that the switch that generates them happens
function probs(star::Int64,fin::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64)
    rates = [ k1, K1*star, k2*(star^3), K2*B*(star^2) ]
    # not actually possible to distinguish two processes, which will effect the probabilities
    if fin - star == 1
        P = (rates[1] + rates[4])/sum(rates)
    elseif fin - star == -1
        P = (rates[2] + rates[3])/sum(rates)
    else
        error()
    end
    return(P)
end

# This function now attempts to differentiate path taken
function probsR(star::Int64,fin::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64)
    rates = [ k1, K1*star, k2*(star^3), K2*B*(star^2) ]
    # not actually possible to distinguish two processes, which will effect the probabilities
    R = rand()
    if fin - star == 1
        if R < rates[1]/(rates[1] + rates[4])
            P = rates[1]/sum(rates)
            rev = 2
        else
            P = rates[4]/sum(rates)
            rev = 3
        end
    elseif fin - star == -1
        if R < rates[2]/(rates[2] + rates[3])
            P = rates[2]/sum(rates)
            rev = 1
        else
            P = rates[3]/sum(rates)
            rev = 4
        end
    else
        error()
    end
    # find reverse rates from final point
    rsR = [ k1, K1*star, k2*(fin^3), K2*B*(fin^2) ]
    Pr = rsR[rev]/sum(rsR)
    return(P,Pr)
end

# This function now just selects most probable of the two
function probsR2(star::Int64,fin::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64)
    rates = [ k1, K1*star, k2*(star^3), K2*B*(star^2) ]
    # not actually possible to distinguish two processes, which will effect the probabilities
    if fin - star == 1
        if rates[1] > rates[4]
            P = rates[1]/sum(rates)
            rev = 2
        else
            P = rates[4]/sum(rates)
            rev = 3
        end
    elseif fin - star == -1
        if rates[2] > rates[3]
            P = rates[2]/sum(rates)
            rev = 1
        else
            P = rates[3]/sum(rates)
            rev = 4
        end
    else
        error()
    end
    # find reverse rates from final point
    rsR = [ k1, K1*star, k2*(fin^3), K2*B*(fin^2) ]
    Pr = rsR[rev]/sum(rsR)
    return(P,Pr)
end

# This function now just selects less probable of the two
function probsR3(star::Int64,fin::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64)
    rates = [ k1, K1*star, k2*(star^3), K2*B*(star^2) ]
    # not actually possible to distinguish two processes, which will effect the probabilities
    if fin - star == 1
        if rates[1] < rates[4]
            P = rates[1]/sum(rates)
            rev = 2
        else
            P = rates[4]/sum(rates)
            rev = 3
        end
    elseif fin - star == -1
        if rates[2] < rates[3]
            P = rates[2]/sum(rates)
            rev = 1
        else
            P = rates[3]/sum(rates)
            rev = 4
        end
    else
        error()
    end
    # find reverse rates from final point
    rsR = [ k1, K1*star, k2*(fin^3), K2*B*(fin^2) ]
    Pr = rsR[rev]/sum(rsR)
    return(P,Pr)
end

function main()
    # first read in langevin trajectories
    one = false
    if one == true
        input_file1 = "../Results/1809/$(ARGS[1])S1.csv"
    else
        input_file1 = "../Results/1809/$(ARGS[1])S2.csv"
    end
    input_filep = "../Results/1809/$(ARGS[1])Sp.csv"
    points1 = Array{Float64,1}(undef,0)
    ps = Array{Float64,1}(undef,0)
    open(input_filep, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            p = parse(Float64, line)
            ps = vcat(ps, p)
        end
    end
    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            X = parse(Float64,line)
            points1 = vcat(points1, X)
        end
    end
    # the parameters will need volume rescaling
    k1 = ps[1]
    K1 = ps[2]
    k2 = ps[3]
    K2 = ps[4]
    B = ps[5]
    Ωi = ps[6] # This is effictively meaningless though
    T1 = ps[7]
    T2 = ps[8]
    Ω = 200000 # new volume
    # rescale rates appropriately
    k1 = k1*Ω
    k2 = k2/(Ω^2)
    K2 = K2/Ω
    # rescale points1
    points1 = points1*Ω
    # Should now make path in A,B
    path = Array{Int64,1}(undef,0)
    ts = Array{Float64,1}(undef,0)
    # approximate steady state
    val = round(Int64,points1[1])
    # initial time t=0
    t = 0
    ts = vcat(ts,t)
    # and add as start of path
    path = vcat(path,val)
    L = length(points1)
    # want to find point where X crossed the 0.5 point
    for i = 2:L
        # first establish if there has been any change from the prior step
        if round(Int64,points1[i]) == path[end]
            # if no change between successive points no reason to update path
        else
            # count changes to A and B between two points
            dX = round(Int64,points1[i]) - path[end]
            th = (i-1)*T1/(L-1) # Time here
            tp = (i-2)*T1/(L-1) # Time at prior step
            vals = Array{Int64,1}(undef,0) # array to store values
            tempt = Array{Float64,1}(undef,0) # vector to temporaily store times
            if dX != 0
                # find difference in X in interval
                dx = (points1[i] - points1[i-1])
                # find first & potentially only point of change
                xp = round(Int64,points1[i]) - 0.5*sign(dX)
                # find distance from point of point of change
                deltaX = points1[i] - xp
                # as a fraction of spacing
                fracX = deltaX/dx
                # find time and save
                t = th - fracX*(th - tp)
                tempt = vcat(tempt,t)
                val = round(Int64,points1[i])
                vals = vcat(vals,val)
                n = 1 # first point done
                # loop until all points are done
                done = false
                while done == false
                    if n == abs(dX)
                        done = true
                        break
                    end
                    deltaX = 1*sign(dX)
                    fracX = deltaX/dx
                    # remove time of gap from time
                    tn = fracX*(th - tp)
                    t -= tn
                    tempt = vcat(tempt,t)
                    val = round(Int64,points1[i]) - n*sign(dX)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            # Now reorder vals by the times
            p = sortperm(tempt)
            tempt = sort(tempt)
            vals = vals[p]
            # then vcat to path
            ts = vcat(ts,tempt)
            path = vcat(path,vals)
        end
    end
    # now need to think about how I use this new path
    # need both forward and backward trajectories in order to establish entropy production
    backp = path[end:-1:1]
    backts = ts[end:-1:1]
    len = length(ts)-1
    Pf = zeros(len)
    Pb = zeros(len)
    for i = 2:len+1
        # Pf[i-1], Pb[i-1] = probsR(path[i-1],path[i],k1,K1,k2,K2,B)
        Pf[i-1] = probs(path[i-1],path[i],k1,K1,k2,K2,B)
        Pb[i-1] = probs(backp[i-1],backp[i],k1,K1,k2,K2,B)
    end
    ents = log.(Pf./Pb[end:-1:1])
    # ents = log.(Pf./Pb)
    plot(ents)
    savefig("../Results/test1.png")
    # average the entropy production
    red = 2500#250
    len2 = floor(Int64,len/red) + 1
    entsr = zeros(len2)
    for i = 1:len2-1
        entsr[i] = sum(ents[(1+(red*(i-1))):(red*i)])
    end
    entsr[end] = sum(ents[(red*(len2-1)+1):end])
    # and make rate of entropy production
    entsp = zeros(len2)
    for i = 1:len2-1
        entsp[i] = entsr[i]/(ts[1+red*i] - ts[1+(red*(i-1))])
    end
    entsp[end] = entsr[end]/(ts[end] - ts[1+(red*(len2-1))])
    # should plot against X here
    X = path[1:red:end,1]/Ω
    plot(X,entsp/Ω,label="Entropy Production")
    scatter!([X[end]],[0],seriescolor=:red,label="")
    scatter!([X[1]],[0],seriescolor=:green,label="")
    plot!(title =  "$(sum(ents)/Ω)")
    if one == true
        savefig("../Results/MastEntsfor.png")
    else
        savefig("../Results/MastEntsbac.png")
    end
    println(sum(ents)/Ω)
    return(nothing)
end

@time main()
