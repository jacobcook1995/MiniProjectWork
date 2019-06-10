# !/usr/bin/env julia
# plotDandS.jl
# A script to read in the data generated by the other scripts and then generate
# appropriate scatter plots for diffusion matrices and steady state entropies
#
# Author: Jacob Cook
# Date: April 2019

# Edit below here
using Plots
using LaTeXStrings
import PyPlot
using Statistics
using Plots.PlotMeasures
using LsqFit

function main()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide a name for toggle switch inputs")
        return(nothing)
    elseif length(ARGS) == 1
        println("Error: Need to provide a name for Schlögl inputs.")
        return(nothing)
    end
    # Check there is a file of Schnakenberg entropy productions to be read
    infile = "../Results/Fig3Data/$(ARGS[1])schnak.csv"
    if ~isfile(infile)
        println("Error: No file of Entropy Productions to be read.")
        return(nothing)
    end
    # now read in Schnakenberg entropy productions
    l = countlines(infile)
    w = 3
    ents = zeros(l,w)
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
                ents[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of Schnakenberg entropy productions to be read for Schlögl
    infile = "../Results/Fig3DataS/$(ARGS[2])schnakS.csv"
    if ~isfile(infile)
        println("Error: No file of Entropy Productions to be read.")
        return(nothing)
    end
    # now read in Schnakenberg entropy productions
    l = countlines(infile)
    w = 3
    Sents = zeros(l,w)
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
                Sents[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now read in toggle switch diffusion matrices
    infile = "../Results/Fig3Data/$(ARGS[1])D.csv"
    if ~isfile(infile)
        println("Error: No file of Diffusion matrices to be read.")
        return(nothing)
    end
    l = countlines(infile)
    w = 3
    DT = zeros(l,w)
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
                DT[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now read in Schlögl diffusion matrices
    infile = "../Results/Fig3DataS/$(ARGS[2])DS.csv"
    if ~isfile(infile)
        println("Error: No file of Diffusion matrices to be read.")
        return(nothing)
    end
    l = countlines(infile)
    w = 3
    DS = zeros(l,w)
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
                DS[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now read in toggle switch entropies
    infile = "../Results/Fig3Data/$(ARGS[1])ent.csv"
    if ~isfile(infile)
        println("Error: No file of Entropies to be read.")
        return(nothing)
    end
    l = countlines(infile)
    w = 2
    entT = zeros(l,w)
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
                entT[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now read in Schlögl entropies
    infile = "../Results/Fig3DataS/$(ARGS[2])entS.csv"
    if ~isfile(infile)
        println("Error: No file of Entropies to be read.")
        return(nothing)
    end
    l = countlines(infile)
    w = 2
    entS = zeros(l,w)
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
                entS[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Line that sets up pyplot and basic sizings
    pyplot(dpi=300,titlefontsize=20,guidefontsize=17,legendfontsize=15,tickfontsize=14)
    # Plot entropy production differences vs diffusion matrix ratio
    scatter([ents[:,3].-ents[:,1]],[DT[:,3]./DT[:,1]],label="",color=1,top_margin=8.0mm)
    plot!(xlabel=L"\dot{S}_{B}-\dot{S}_{A}",ylabel=L"mag(D_{B})/mag(D_{A})",title=L"mag(D_{B})/mag(D_{A})\;vs\;\Delta\dot{S}")
    scatter!([Sents[:,3].-Sents[:,1]],[DS[:,3]./DS[:,1]],label="",color=2)
    savefig("../Results/DiffDvsEntProd.png")
    # plot steady state entropies against entropy productions
    scatter([ents[:,1].-ents[:,3]],[entT[:,1].-entT[:,2]],label="",ylim=(-10,10),color=1)
    scatter!([Sents[:,1].-Sents[:,3]],[entS[:,1].-entS[:,2]],label="",color=2)
    plot!(xlabel=L"\dot{S}_A-\dot{S}_B",ylabel=L"S_A - S_B",title=L"\Delta S\;vs\;\Delta\dot{S}",top_margin=8.0mm)
    savefig("../Results/DiffEnt.png")
    return(nothing)
end

@time main()