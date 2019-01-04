#!/usr/bin/env julia
# EntProdtraj.jl
# A script to read in parameters and trajcetories and then to generate nice plots
# of the entropy productions by two differenet methods.
# 1) using the langevin entropy production term
# 2) by using a large volume reduced form of the master equation
#
# Author: Jacob Cook
# Date: January 2019

using Plots
using LaTeXStrings
using PyCall
import PyPlot

# Diffusion matrix D
# ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
function D!(D::Array{Float64,2},x::Array{Float64,1},ps::Array{Float64,1})
    D[1,1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) + (ps[5]+ps[2])*x[1] + ps[6]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) + (ps[7]+ps[4])*x[2] + ps[8]
    return(D)
end

# vector of forces
function f!(f::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) - ps[2]*x[1] - ps[5]*x[1] + ps[6]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) - ps[4]*x[2] - ps[7]*x[2] + ps[8]
    return(f)
end

# A function to find and plot the langevin entropy productions of the various trajectories
function LangEnt(traj::Array{Float64,2},ps::Array{Float64,1},dt::Float64)
    entp = zeros(size(traj,1)-1)
    f = [0.0; 0.0]
    D = [0.0 0.0; 0.0 0.0]
    for i = 1:length(entp)
        qdot = (traj[i+1,:] .- traj[i,:])/(dt)
        # and need to find midpoint
        posA = (traj[i+1,1] + traj[i,1])/(2)
        posB = (traj[i+1,2] + traj[i,2])/(2)
        f = f!(f,[posA,posB],ps)
        D = D!(D,[posA,posB],ps)
        for j = 1:2
            for k = 1:2
                if D[j,k] != 0
                    entp[i] += 2*qdot[j]*f[k]*dt/D[j,k]
                end
            end
        end
    end
    return(entp)
end

# function that takes in two points and finds the probability that the switch that generates them happens
function probs(star::Array{Int64,1},fin::Array{Int64,1},k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ k*r/(r+f*star[2]^2), kmin*star[1], K*star[1], Kmin, q*r/(r+f*star[1]^2), qmin*star[2], Q*star[2], Qmin ]
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

# A function to find and plot the reduced master equation entropy productions of the various trajectories
function MastEnt(traj::Array{Float64,2},ps::Array{Float64,1},T::Float64,Ω::Int64,red::Int64)
    k = ps[1]
    kmin = ps[2]
    q = ps[3]
    qmin = ps[4]
    K = ps[5]
    Kmin = ps[6]
    Q = ps[7]
    Qmin = ps[8]
    r = ps[9]
    f = ps[10]
    # rescale rates appropriately
    k = k*Ω
    Kmin = Kmin*Ω
    q = q*Ω
    Qmin = Qmin*Ω
    f = f/(Ω^2)
    # rescale trajectory
    traj = traj*Ω
    # Should now make path in A,B
    path = Array{Int64,2}(undef,0,2)
    ts = Array{Float64,1}(undef,0)
    # approximate steady state
    val = fill(0,1,2)
    for i = 1:2
        val[i] = round(Int64,traj[1,i])
    end
    # initial time t=0
    t = 0
    ts = vcat(ts,t)
    # and add as start of path
    path = vcat(path,val)
    L = size(traj,1)
    # want to find point where A/B crossed the 0.5 point
    for i = 2:L
        # first establish if there has been any change from the prior step
        if round(Int64,traj[i,1]) == path[end,1] && round(Int64,traj[i,2]) == path[end,2]
            # if no change between successive points no reason to update path
        else
            # count changes to A and B between two points
            dA = round(Int64,traj[i,1]) - path[end,1]
            dB = round(Int64,traj[i,2]) - path[end,2]
            th = (i-1)*T/(L-1) # Time here
            tp = (i-2)*T/(L-1) # Time at prior step
            vals = Array{Int64,2}(undef,0,2) # array to store values
            tempt = Array{Float64,1}(undef,0) # vector to temporaily store times
            if dA != 0
                # find difference in A in interval
                da = (traj[i,1] - traj[i-1,1])
                # find first & potentially only point of change
                ap = round(Int64,traj[i,1]) - 0.5*sign(dA)
                # find distance from point of point of change
                deltaA = traj[i,1] - ap
                # as a fraction of spacing
                fracA = deltaA/da
                # find time and save
                t = th - fracA*(th - tp)
                tempt = vcat(tempt,t)
                # find B and save point
                B = traj[i,2] - fracA*(traj[i,2] - traj[i-1,2])
                val[1] = round(Int64,traj[i,1])
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
                    B -= fracA*(traj[i,2] - traj[i-1,2])
                    val[1] = round(Int64,traj[i,1]) - n*sign(dA)
                    val[2] = round(Int64,B)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            if dB != 0
                # find difference in B in interval
                db = (traj[i,2] - traj[i-1,2])
                # find first & potentially only point of change
                bp = round(Int64,traj[i,2]) - 0.5*sign(dB)
                # find distance from point of point of change
                deltaB = traj[i,2] - bp
                # as a fraction of spacing
                fracB = deltaB/db
                # find time and save
                t = th - fracB*(th - tp)
                tempt = vcat(tempt,t)
                # find A and save point
                A = traj[i,1] - fracB*(traj[i,1] - traj[i-1,1])
                val[1] = round(Int64,A)
                val[2] = round(Int64,traj[i,2])
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
                    A -= fracB*(traj[i,1] - traj[i-1,1])
                    val[1] = round(Int64,A)
                    val[2] = round(Int64,traj[i,2]) - n*sign(dB)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            # Now reorder vals by the times
            p = sortperm(tempt)
            tempt = sort(tempt)
            vals = vals[p,:]
            # then vcat to path
            ts = vcat(ts,tempt)
            path = vcat(path,vals)
        end
    end
    # need both forward and backward trajectories in order to establish entropy production
    backp = path[end:-1:1,:]
    backts = ts[end:-1:1]
    len = length(ts)-1
    Pf = zeros(len)
    Pb = zeros(len)
    for i = 2:len+1
        Pf[i-1] = probs(path[i-1,:],path[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        Pb[i-1] = probs(backp[i-1,:],backp[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
    end
    entp = log.(Pf./Pb[end:-1:1])
    len2 = floor(Int64,len/red)+1
    entsr = zeros(len2)
    As = zeros(len2)
    for i = 1:len2-1
        As[i] = sum(path[(1+(red*(i-1))):(red*i),1])/red
        entsr[i] = sum(entp[(1+(red*(i-1))):(red*i)])
    end
    entsr[end] = sum(entp[(red*(len2-1)+1):end])
    As[end] = sum(path[(red*(len2-1)+1):end,1])/length(path[(red*(len2-1)+1):end,1])
    # finally extract the path in A
    entp = entsr
    return(entp,As)
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
    infile = "../Results/Fig2Data/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    len = countlines(infile)
    w = 10
    ps = zeros(len,w)
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
    # now read in the 6 trajcetories
    N = 600
    traj = zeros(N+1,2*2*len)
    for i = 1:len
        for j = 1:2
            if j == 1
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])A2B.csv"
            else
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])B2A.csv"
            end
            w = 2
            open(infile, "r") do in_file
                # Use a for loop to process the rows in the input file one-by-one
                k = 1
                for line in eachline(in_file)
                    # parse line by finding commas
                    L = length(line)
                    comma = fill(0,w+1)
                    m = 1
                    for l = 1:L
                        if line[l] == ','
                            m += 1
                            comma[m] = l
                        end
                    end
                    comma[end] = L+1
                    for l = 1:w
                        traj[k,4*(i-1)+2*(j-1)+l] = parse(Float64,line[(comma[l]+1):(comma[l+1]-1)])
                    end
                    k += 1
                end
            end
        end
    end
    # Now read in additional data from the other stuff
    w = 4
    Act = zeros(2*len,4)
    for i = 1:len
        for j = 1:2
            if j == 1
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])A2BD.csv"
            else
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])B2AD.csv"
            end
            open(infile, "r") do in_file
                # Single line file so doesn't matter
                for line in eachline(in_file)
                    # parse line by finding commas
                    L = length(line)
                    comma = fill(0,w+1)
                    m = 1
                    for l = 1:L
                        if line[l] == ','
                            m += 1
                            comma[m] = l
                        end
                    end
                    comma[end] = L+1
                    for l = 1:w
                        Act[2*(i-1)+j,l] = parse(Float64,line[(comma[l]+1):(comma[l+1]-1)])
                    end
                end
            end
        end
    end
    # Check there is a file of steady states to be read
    infile = "../Results/Fig2Data/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    len = countlines(infile)
    w = 6
    steads = zeros(len,w)
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
    # Now calculate and plot langevin entropy productions
    LatS = L"\Delta S_{L}"
    pyplot()
    # for i = 1:len
    #     for j = 1:2
    #         d = 2*(j-1)+4*(i-1)
    #         entp = LangEnt(traj[:,(1+d):(2+d)],ps[i,:],Act[2*(i-1)+j,1]/N)
    #         plot(traj[1:end-1,1+d],entp,label=LatS,dpi=300,legend=:best)
    #         plot!(xlabel="A",ylabel="Ent Prod")
    #         scatter!([steads[i,1+4*(j-1)]],[0.0],color=:green,label="Start")
    #         scatter!([steads[i,3]],[0.0],color=:yellow,label="Saddle")
    #         scatter!([steads[i,5-4*(j-1)]],[0.0],color=:red,label="End")
    #         savefig("../Results/Fig2Graphs/$(i)$(j)LangEnt.png")
    #     end
    # end
    LatS2 = L"\Delta S_{M}"
    Ω = 5000
    LatS3 = latexstring("\\Omega = $(Ω)")
    red = 250 # batch size
    for i = 1:len
        for j = 1:2
            d = 2*(j-1)+4*(i-1)
            entp2, As = MastEnt(traj[:,(1+d):(2+d)],ps[i,:],Act[2*(i-1)+j,1],Ω,red)
            plot(As/Ω,entp2/Ω,label=LatS2,dpi=300,legend=:best)
            plot!(xlabel="A",ylabel="Ent Prod",legendtitle=LatS3)
            scatter!([steads[i,1+4*(j-1)]],[0.0],color=:green,label="Start")
            scatter!([steads[i,3]],[0.0],color=:yellow,label="Saddle")
            scatter!([steads[i,5-4*(j-1)]],[0.0],color=:red,label="End")
            savefig("../Results/Fig2Graphs/$(i)$(j)MastEnt.png")
        end
    end
    return(nothing)
end

@time main()
