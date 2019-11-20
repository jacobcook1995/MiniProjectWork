# !/usr/bin/env julia
# epef.jl
# A script that that makes graphs comparing the entropy production (EF) and entropy
# flow (EF) terms of both models. A large number of plots shall be generated for the comparisons
#
# Author: Jacob Cook
# Date: November 2019

using Plots
using LaTeXStrings
using Statistics
import PyPlot

# One plotting function in this script that should hopefully obtain everything needed
function first()
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
    # Now ready to do first load of plots
    pyplot(dpi=200)#,titlefontsize=17,guidefontsize=14,legendfontsize=15,tickfontsize=14)
    scatter(ent[:,1:2],ent[:,3:4],label=[L"\dot{S}_A" L"\dot{S}_B"])
    plot!(xlabel="Exact (Schnakenberg)",ylabel="EP term",title="Steady state entropy productions toggle")
    savefig("../Results/EPEF/steadT.png")
    scatter(Sent[:,1:2],Sent[:,3:4],label=[L"\dot{S}_A" L"\dot{S}_B"])
    plot!(xlabel="Exact (Schnakenberg)",ylabel="EP term",title="Steady state entropy productions Schlögl")
    savefig("../Results/EPEF/steadS.png")
    scatter(ent[:,1].-ent[:,2],ent[:,3].-ent[:,4],label=L"\dot{S}_A - \dot{S}_B")
    plot!(xlabel="Exact (Schnakenberg)",ylabel="EP term",title="Diff. Steady state entropy productions toggle")
    savefig("../Results/EPEF/diffsteadT.png")
    scatter(Sent[:,1].-Sent[:,2],Sent[:,3].-Sent[:,4],label=L"\dot{S}_A - \dot{S}_B")
    plot!(xlabel="Exact (Schnakenberg)",ylabel="EP term",title="Diff. Steady state entropy productions Schlögl")
    savefig("../Results/EPEF/diffsteadS.png")
    return(nothing)
end

# A second function to find production and flow terms along the paths
function second()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide a name for toggle switch inputs")
        return(nothing)
    elseif length(ARGS) == 1
        println("Error: Need to provide a name for Schlögl inputs.")
        return(nothing)
    end
    # Can find the EP and EF values from a file
    # Make a structure to store the EP and EF terms for the toggle switch model
    l = 100
    ef = zeros(l,8)
    datayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2Bpf.csv"
            else
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2Apf.csv"
            end
            # check it exits
            if isfile(infile)
                w = 4
                open(infile, "r") do in_file
                    # only one line now
                    for line in eachline(in_file)
                        # parse line by finding commas
                        L = length(line)
                        comma = fill(0,w+1)
                        m = 1
                        for k = 1:L
                            if line[k] == ','
                                m += 1
                                comma[m] = k
                            end
                        end
                        comma[end] = L+1
                        for k = 1:w
                            ef[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
                        end
                    end
                end
            else # if file doesn't exist inform user of missing data and exclude from plots
                if j == 1
                    println("Missing data for run $(i) high A to high B")
                    datayn[i] = false
                else
                    println("Missing data for run $(i) high B to high A")
                    datayn[i] = false
                end
            end
        end
    end
    # Make a structure to store the EP and EF terms for the Schlögl model
    Sef = zeros(l,8)
    Sdatayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])h2lpf.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])l2hpf.csv"
            end
            # check it exits
            if isfile(infile)
                w = 4
                open(infile, "r") do in_file
                    # only one line now
                    for line in eachline(in_file)
                        # parse line by finding commas
                        L = length(line)
                        comma = fill(0,w+1)
                        m = 1
                        for k = 1:L
                            if line[k] == ','
                                m += 1
                                comma[m] = k
                            end
                        end
                        comma[end] = L+1
                        for k = 1:w
                            Sef[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
                        end
                    end
                end
            else # if file doesn't exist inform user of missing data and exclude from plots
                if j == 1
                    println("Missing data for run $(i) high A to high B")
                    Sdatayn[i] = false
                else
                    println("Missing data for run $(i) high B to high A")
                    Sdatayn[i] = false
                end
            end
        end
    end
    # Now need to obtain the master equation switching data
    means = zeros(l,4)
    for i = 1:l
        infile1 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2BMast.csv"
        infile2 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2AMast.csv"
        # Extract forward data
        ne = 101
        w = ne
        mastf = zeros(w-1)
        open(infile1, "r") do in_file
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
                for i = 1:w-1
                    mastf[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Then extract backwards data
        mastb = zeros(w-1)
        open(infile2, "r") do in_file
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
                for i = 1:w-1
                    mastb[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Now need to extract averages and standard deviations and plot against other data
        means[i,1] = mean(mastf)
        means[i,2] = mean(mastb)
        means[i,3] = stdm(mastf,means[i,1];corrected=true)/sqrt(ne-1)
        means[i,4] = stdm(mastb,means[i,2];corrected=true)/sqrt(ne-1)
    end
    # Now repeat entire process for the Schlögl model
    Smeans = zeros(l,4)
    for i = 1:l
        infile1 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])h2lMast.csv"
        infile2 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])l2hMast.csv"
        # Extract forward data
        ne = 101
        w = ne
        mastf = zeros(w-1)
        open(infile1, "r") do in_file
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
                for i = 1:w-1
                    mastf[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Then extract backwards data
        mastb = zeros(w-1)
        open(infile2, "r") do in_file
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
                for i = 1:w-1
                    mastb[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Now need to extract averages and standard deviations and plot against other data
        Smeans[i,1] = mean(mastf)
        Smeans[i,2] = mean(mastb)
        Smeans[i,3] = stdm(mastf,means[i,1];corrected=true)/sqrt(ne-1)
        Smeans[i,4] = stdm(mastb,means[i,2];corrected=true)/sqrt(ne-1)
    end
    # First compare EP terms with master equation entropy productions
    pyplot(dpi=200)
    plot(xlabel="EP term",ylabel="Exact (Gillespie)",title="Path entropy productions toggle")
    # First plot one set
    for i = 1:l
        if i == 1
            plot!([ef[i,1]],[means[i,1]],yerror=[means[i,3]],label=L"\Delta S_{A{\rightarrow}B}",color=1)
        else
            plot!([ef[i,1]],[means[i,1]],yerror=[means[i,3]],label="",color=1)
        end
    end
    # and then the other
    for i = 1:l
        if i == 1
            plot!([ef[i,5]],[means[i,2]],yerror=[means[i,4]],label=L"\Delta S_{B{\rightarrow}A}",color=2)
        else
            plot!([ef[i,5]],[means[i,2]],yerror=[means[i,4]],label="",color=2)
        end
    end
    savefig("../Results/EPEF/switchT.png")
    plot(xlabel="EP term",ylabel="Exact (Gillespie)",title="Path entropy productions Schlögl")
    # First plot one set
    for i = 1:l
        if i == 1
            plot!([Sef[i,1]],[Smeans[i,1]],yerror=[Smeans[i,3]],label=L"\Delta S_{A{\rightarrow}B}",color=1)
        else
            plot!([Sef[i,1]],[Smeans[i,1]],yerror=[Smeans[i,3]],label="",color=1)
        end
    end
    # and then the other
    for i = 1:l
        if i == 1
            plot!([Sef[i,5]],[Smeans[i,2]],yerror=[Smeans[i,4]],label=L"\Delta S_{B{\rightarrow}A}",color=2)
        else
            plot!([Sef[i,5]],[Smeans[i,2]],yerror=[Smeans[i,4]],label="",color=2)
        end
    end
    savefig("../Results/EPEF/switchS.png")
    # Do differences
    # Then same with entropy flow
    # Also need to compare to Langevin entropy production => As Robert is sure to ask about this
    return(nothing)
end

@time second()
