# !/usr/bin/env julia
# epef.jl
# A script that that makes graphs comparing the entropy production (EF) and entropy
# flow (EF) terms of both models. A large number of plots shall be generated for the comparisons
#
# Author: Jacob Cook
# Date: November 2019

using Plots
using LaTeXStrings
import PyPlot

# One main plotting function in this script that should hopefully obtain everything needed
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
    # Firstly wish to make steady state comparisons
    # read in EP terms => as EF terms equal and opposite at steady state
    # Check there is a file of productions to be read
    infile = "../Results/Fig3Data/$(ARGS[1])prod.csv"
    if ~isfile(infile)
        println("Error: No file of EP terms to be read for toggle switch.")
        return(nothing)
    end
    # now do readings
    l = countlines(infile)
    w = 3
    prods = zeros(l,w)
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
                prods[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of EP's to be read for Schlögl model
    infile = "../Results/Fig3DataS/$(ARGS[2])prodS.csv"
    if ~isfile(infile)
        println("Error: No file of EP terms to be read for Schlögl model.")
        return(nothing)
    end
    # now read in EP terms for
    l = countlines(infile)
    w = 3
    Sprods = zeros(l,w)
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
                Sprods[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of steady-state entropy productions to be read
    infile = "../Results/Fig3Data/$(ARGS[1])schnak.csv"
    if ~isfile(infile)
        println("Error: No file of Entropy Productions to be read for toggle switch.")
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
    # Check there is a file of steady-state entropy productions to be read for Schlögl
    infile = "../Results/Fig3DataS/$(ARGS[2])schnakS.csv"
    if ~isfile(infile)
        println("Error: No file of Entropy Productions to be read for Schlögl.")
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
    # combine entropy production and EP terms into one structure for toggle
    ent = zeros(l,4)
    for i = 1:l
        ent[i,1] = prods[i,1] # A
        ent[i,2] = prods[i,3] # B
        ent[i,3] = ents[i,1] # A
        ent[i,4] = ents[i,3] # B
    end
    # combine entropy production and EP into one structure for Schlögl model
    Sent = zeros(l,4)
    for i = 1:l
        Sent[i,1] = Sprods[i,1] # A
        Sent[i,2] = Sprods[i,3] # B
        Sent[i,3] = Sents[i,1] # A
        Sent[i,4] = Sents[i,3] # B
    end
    return(nothing)
end

@time main()
