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
    for i = 1:len
        for j = 1:2
            d = 2*(j-1)+4*(i-1)
            entp = LangEnt(traj[:,(1+d):(2+d)],ps[i,:],Act[2*(i-1)+j,1]/N)
            plot(traj[1:end-1,1+d],entp,label=LatS,dpi=300,legend=:best)
            plot!(xlabel="A",ylabel="Ent Prod")
            scatter!([steads[i,1+4*(j-1)]],[0.0],color=:green,label="Start")
            scatter!([steads[i,3]],[0.0],color=:yellow,label="Saddle")
            scatter!([steads[i,5-4*(j-1)]],[0.0],color=:red,label="End")
            savefig("../Results/Fig2Graphs/$(i)$(j)LangEnt.png")
        end
    end
    LatS2 = L"\\Delta S_{M}"
    return(nothing)
end

@time main()
