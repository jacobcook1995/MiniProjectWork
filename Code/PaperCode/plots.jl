#!/usr/bin/env julia
# plots.jl
# A script to read in the data generated by the other scripts and then generate
# appropriate scatter plots
#
# Author: Jacob Cook
# Date: November 2018

using Plots
using LaTeXStrings
using PyCall
pygui(:qt5)
import PyPlot

function main()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])prod.csv"
    if ~isfile(infile)
        println("Error: No file of 'Entropy Productions' to be read.")
        return(nothing)
    end
    # now read in 'Entropy productions'
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
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])gill.csv"
    if ~isfile(infile)
        println("Error: No file of 'Entropy Productions' to be read.")
        return(nothing)
    end
    # now read in 'Entropy productions'
    l = countlines(infile)
    w = 2
    tjents = zeros(l,w)
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
                tjents[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])schnak.csv"
    if ~isfile(infile)
        println("Error: No file of Entropy Productions to be read.")
        return(nothing)
    end
    # now read in parameters
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
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])ent.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 2
    shent = zeros(l,w)
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
                shent[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])D.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 3
    Ds = zeros(l,w)
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
                Ds[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # combine entropy productions as one structure
    ent = zeros(l,6)
    for i = 1:l
        ent[i,1] = prods[i,1]
        ent[i,2] = prods[i,3]
        ent[i,3] = ents[i,1]
        ent[i,4] = ents[i,3]
        ent[i,5] = tjents[i,1]
        ent[i,6] = tjents[i,2]
    end
    # find points where
    wrong = Array{Int64,1}(undef,0)
    for i = 1:l
        if ent[i,1] < ent[i,2] && ent[i,5] > ent[i,6]
            wrong = vcat(wrong,i)
        elseif ent[i,1] > ent[i,2] && ent[i,5] < ent[i,6]
            wrong = vcat(wrong,i)
        end
    end
    pyplot()
    scatter([ent[:,3]],[Ds[:,1]],label="")
    scatter!([ent[:,4]],[Ds[:,3]],label="")
    plot!(xlabel=L"\dot{S}",ylabel=L"mag(D)",title="Magnitude of D vs Entropy Production")
    savefig("../Results/DvsEntProd.png")
    scatter([ent[:,3].-ent[:,4]],[Ds[:,1]./Ds[:,3]],label="")
    plot!(xlabel=L"\dot{S}_{h}-\dot{S}_{l}",ylabel=L"mag(D_{h})/mag(D_{l})",title="Ratio of D vs Entropy Production diff")
    savefig("../Results/DiffDvsEntProd.png")
    # scatter([ent[:,1].-ent[:,2]],[ent[:,3].-ent[:,4]],label="")
    # plot!(xlabel=L"\dot{S}_{h,prod}-\dot{S}_{l,prod}",ylabel=L"\dot{S}_h - \dot{S}_l",title="Diff in Entropy vs Production")
    # savefig("../Results/DifferProds.png")
    # # Plot entropies against Schnakenberg entropy production
    # scatter([shent[:,1]],[ent[:,3]],label="")
    # scatter!([shent[:,2]],[ent[:,4]],label="")
    # plot!(xlabel=L"\dot{S}",ylabel=L"S",title="Ent vs Ent Prod")
    # savefig("../Results/EntvsSchnak.png")
    # scatter([shent[:,1].-shent[:,2]],[ent[:,3].-ent[:,4]],label="")
    # plot!(xlabel=L"\dot{S}_h-\dot{S}_l",ylabel=L"S_h - S_l",title="Diff in Entropy vs Production")
    # savefig("../Results/DiffEnt.png")
    # plot scatter graph of the wrong points
    # plot(xlabel="Ent Prod Rate Term",ylabel="Ent Prod Rate Traj",title="Mismatched points")
    # for i = 1:length(wrong)
    #     scatter!([ent[wrong[i],1],ent[wrong[i],2]], [ent[wrong[i],5],ent[wrong[i],6]],label="")
    # end
    # savefig("../Results/WrongTrajvsTerms.png")
    # then plot scatter graph
    # scatter([ent[:,1]], [ent[:,3]],label="")
    # scatter!([ent[:,2]], [ent[:,4]],label="")
    # plot!(xlabel="Ent Prod Rate Term",ylabel="Ent Prod Rate Sch",title="Ent Prod Terms vs Schnakenberg")
    # savefig("../Results/SchnakvsTerms.png")
    # scatter([ent[:,1]], [ent[:,5]],label="")
    # scatter!([ent[:,2]], [ent[:,6]],label="")
    # plot!(xlabel="Ent Prod Rate Term",ylabel="Ent Prod Rate Traj",title="Trajectory Entropy Prod vs Ent Prod Terms")
    # savefig("../Results/TermsvsTraj.png")
    # scatter([ent[:,3]], [ent[:,5]],label="")
    # scatter!([ent[:,4]], [ent[:,6]],label="")
    # plot!(xlabel="Ent Prod Rate Sch",ylabel="Ent Prod Rate Traj",title="Trajectory Entropy Prod vs Schnakenberg")
    # savefig("../Results/SchnakvsTraj.png")
    return(nothing)
end

@time main()
