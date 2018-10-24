#!/usr/bin/env julia
# trajent.jl
# A script to take in a trajectory from the Langevin formulalism convert to a sequence
# of jumps in the master equation formulation, and then calculate the entropy production of this process
#
# Author: Jacob Cook
# Date: October 2018
using Plots

function main()
    # first read in langevin trajectories
    input_file1 = "../Results/1809/$(ARGS[1])1.csv"
    input_file2 = "../Results/1809/$(ARGS[1])2.csv"
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
    # now do second file
    open(input_file2, "r") do in_file
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
            points2 = vcat(points2, [ A B ])
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
    Ω = 10000 # new volume
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
    plot(path[:,1],path[:,2])
    savefig("../Results/test1.png")
    plot(ts)
    savefig("../Results/test2.png")
    return(nothing)
end

@time main()
