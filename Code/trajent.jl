#!/usr/bin/env julia
# trajent.jl
# A script to take in a trajectory from the Langevin formulalism convert to a sequence
# of jumps in the master equation formulation, and then calculate the entropy production of this process
#
# Author: Jacob Cook
# Date: October 2018
using Plots

# function that takes in two points and finds the probability that the switch that generates them happens
function probs(star::Array{Int64,1},fin::Array{Int64,1},k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ k*r/(r+f*star[2]^2), kmin*star[1], K*star[1], Kmin, q*r/(r+f*star[1]^2), qmin*star[2], Q*star[2], Qmin]
    # not actually possible to distinguish two processes, which will effect the probabilities
    if fin[1] - star[1] == 0
        if fin[2] - star[2] == 1
            P = (rates[5] + rates[8])/sum(rates)
        elseif fin[2] - star[2] == -1
            P = (rates[6] + rates[7])/sum(rates)
        else
            error()
        end
    else
        if fin[1] - star[1] == 1
            P = (rates[1] + rates[4])/sum(rates)
        elseif fin[1] - star[1] == -1
            P = (rates[2] + rates[3])/sum(rates)
        else
            error()
        end
    end
    return(P)
end

# Also need another function that calculates probability of waiting time t
# function that takes in two points and finds the probability that the switch that generates them happens
function probw(star::Array{Int64,1},fin::Array{Int64,1},k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,t1::Float64,t2::Float64)
    rates = [ k*r/(r+f*star[2]^2), kmin*star[1], K*star[1], Kmin, q*r/(r+f*star[1]^2), qmin*star[2], Q*star[2], Qmin]
    t = convert(BigFloat,abs(t2 - t1)) # time waited
    P1 = exp(-t*sum(rates)) # probability of survining to time t
    ϵ = 10.0^-15
    P2 = exp(-t*sum(rates)) - exp(-(t+ϵ)*sum(rates))
    P = P1*P2
    return(P)
end

function main()
    # first read in langevin trajectories
    one = true
    if one == true
        input_file1 = "../Results/1809/$(ARGS[1])1.csv"
    else
        input_file1 = "../Results/1809/$(ARGS[1])2.csv"
    end
    input_filep = "../Results/1809/$(ARGS[1])p.csv"
    points1 = Array{Float64,2}(undef,0,2)
    points2 = Array{Float64,2}(undef,0,2)
    ps = Array{Float64,1}(undef,0)
    open(input_filep, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            p = parse(Float64, line)
            ps = vcat(ps, p)
        end
    end
    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for j = 1:L
                if line[j] == ','
                    comma = j
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points1 = vcat(points1, [ A B ])
        end
    end
    # the parameters will need volume rescaling
    K = ps[1]
    k = ps[2]
    Q = ps[3]
    q = ps[4]
    kmin = ps[5]
    Kmin = ps[6]
    qmin = ps[7]
    Qmin = ps[8]
    f = ps[9]
    r = ps[10]
    T1 = ps[11]
    T2 = ps[12]
    Ωi = 2 # orginal volume here
    Ω = 25000 # new volume
    # rescale rates appropriately
    k = k*(Ω/Ωi)
    Kmin = Kmin*(Ω/Ωi)
    q = q*(Ω/Ωi)
    Qmin = Qmin*(Ω/Ωi)
    f = f/((Ω/Ωi)^2)
    # rescale points1 and 2
    points1 = points1*(Ω/Ωi)
    points2 = points2*(Ω/Ωi)
    # Should now make path in A,B
    path = Array{Int64,2}(undef,0,2)
    ts = Array{Float64,1}(undef,0)
    # approximate steady state
    val = fill(0,1,2)
    for i = 1:2
        val[i] = round(Int64,points1[1,i])
    end
    # initial time t=0
    t = 0
    ts = vcat(ts,t)
    # and add as start of path
    path = vcat(path,val)
    L = size(points1,1)
    # want to find point where A/B crossed the 0.5 point
    for i = 2:L
        # first establish if there has been any change from the prior step
        if round(Int64,points1[i,1]) == path[end,1] && round(Int64,points1[i,2]) == path[end,2]
            # if no change between successive points no reason to update path
        else
            # count changes to A and B between two points
            dA = round(Int64,points1[i,1]) - path[end,1]
            dB = round(Int64,points1[i,2]) - path[end,2]
            th = (i-1)*T1/(L-1) # Time here
            tp = (i-2)*T1/(L-1) # Time at prior step
            vals = Array{Int64,2}(undef,0,2) # array to store values
            tempt = Array{Float64,1}(undef,0) # vector to temporaily store times
            if dA != 0
                # find difference in A in interval
                da = (points1[i,1] - points1[i-1,1])
                # find first & potentially only point of change
                ap = round(Int64,points1[i,1]) - 0.5*sign(dA)
                # find distance from point of point of change
                deltaA = points1[i,1] - ap
                # as a fraction of spacing
                fracA = deltaA/da
                # find time and save
                t = th - fracA*(th - tp)
                tempt = vcat(tempt,t)
                # find B and save point
                B = points1[i,2] - fracA*(points1[i,2] - points1[i-1,2])
                val[1] = round(Int64,points1[i,1])
                val[2] = round(Int64,B)
                vals = vcat(vals,val)
                n = 1 # first point done
                # loop until all points are done
                done = false
                while done == false
                    if n == abs(dA)
                        done = true
                        break
                    end
                    deltaA = 1*sign(dA)
                    fracA = deltaA/da
                    # remove time of gap from time
                    tn = fracA*(th - tp)
                    t -= tn
                    tempt = vcat(tempt,t)
                    B -= fracA*(points1[i,2] - points1[i-1,2])
                    val[1] = round(Int64,points1[i,1]) - n*sign(dA)
                    val[2] = round(Int64,B)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            if dB != 0
                # find difference in B in interval
                db = (points1[i,2] - points1[i-1,2])
                # find first & potentially only point of change
                bp = round(Int64,points1[i,2]) - 0.5*sign(dB)
                # find distance from point of point of change
                deltaB = points1[i,2] - bp
                # as a fraction of spacing
                fracB = deltaB/db
                # find time and save
                t = th - fracB*(th - tp)
                tempt = vcat(tempt,t)
                # find A and save point
                A = points1[i,1] - fracB*(points1[i,1] - points1[i-1,1])
                val[1] = round(Int64,A)
                val[2] = round(Int64,points1[i,2])
                vals = vcat(vals,val)
                n = 1 # first point done
                # loop until all points are done
                done = false
                while done == false
                    if n == abs(dB)
                        done = true
                        break
                    end
                    deltaB = 1*sign(dB)
                    fracB = deltaB/db
                    # remove time of gap from time
                    tn = fracB*(th - tp)
                    t -= tn
                    tempt = vcat(tempt,t)
                    A -= fracB*(points1[i,1] - points1[i-1,1])
                    val[1] = round(Int64,A)
                    val[2] = round(Int64,points1[i,2]) - n*sign(dB)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            # times are non-sensenscial for cases where A and B both run
            # Now reorder vals by the times
            p = sortperm(tempt)
            tempt = sort(tempt)
            vals = vals[p,:]
            # then vcat to path
            ts = vcat(ts,tempt)
            path = vcat(path,vals)
        end
    end
    # now need to think about how I use this new path
    # need both forward and backward trajectories in order to establish entropy production
    backp = path[end:-1:1,:]
    backts = ts[end:-1:1]
    len = length(ts)-1
    Pf = zeros(len)
    Pb = zeros(len)
    Pwf = zeros(BigFloat,len)
    Pwb = zeros(BigFloat,len)
    Pwf1 = zeros(len)
    Pwb1 = zeros(len)
    Pwf2 = zeros(len)
    Pwb2 = zeros(len)
    for i = 2:len+1
        Pf[i-1] = probs(path[i-1,:],path[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        Pb[i-1] = probs(backp[i-1,:],backp[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        Pwf[i-1] = probw(path[i-1,:],path[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,ts[i-1],ts[i])
        Pwb[i-1] = probw(backp[i-1,:],backp[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,backts[i-1],backts[i])
    end
    ents = log.(Pf./Pb[end:-1:1])
    ents2 = log.(Pwf./Pwb[end:-1:1])
    # plot(ents)
    # savefig("../Results/test1.png")
    # plot(ents2)
    # savefig("../Results/test2.png")
    # average the entropy production
    red = 250
    len2 = floor(Int64,len/red)+1
    entsr = zeros(len2)
    ents2r = zeros(len2)
    for i = 1:len2-1
        entsr[i] = sum(ents[(1+(red*(i-1))):(red*i)])
        ents2r[i] = sum(ents2[(1+(red*(i-1))):(red*i)])
    end
    entsr[end] = sum(ents[(red*(len2-1)+1):end])
    ents2r[end] = sum(ents2[(red*(len2-1)+1):end])
    # and make rate of entropy production
    entsp = zeros(len2)
    for i = 1:len2-1
        entsp[i] = entsr[i]/(ts[1+red*i] - ts[1+(red*(i-1))])
    end
    entsp[end] = entsr[end]/(ts[end] - ts[1+(red*(len2-1))])
    # should plot against A here
    A = path[1:red:end,1]*(Ωi/Ω)
    plot(A,entsp*(Ωi/Ω),label="Entropy Production")
    scatter!([A[end]],[0],seriescolor=:red,label="")
    scatter!([A[1]],[0],seriescolor=:green,label="")
    plot!(title =  "$(sum(ents)*(Ωi/Ω))")
    if one == true
        savefig("../Results/MastEntsB.png")
    else
        savefig("../Results/MastEntsA.png")
    end
    println(sum(ents)*(Ωi/Ω))
    println(sum(log.(Pwf./Pwb[end:-1:1]))*(Ωi/Ω))
    return(nothing)
end

@time main()
