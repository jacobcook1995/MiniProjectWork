# !/usr/bin/env julia
# combplots.jl
# A script to read in the data generated by the other scripts and then generate
# appropriate scatter plots
#
# Author: Jacob Cook
# Date: February 2019

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
    # Check there is a file of productions to be read
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
    # Check there is a file of productions to be read for Schlögl model
    infile = "../Results/Fig3DataS/$(ARGS[2])prodS.csv"
    if ~isfile(infile)
        println("Error: No file of 'Entropy Productions' to be read.")
        return(nothing)
    end
    # now read in 'Entropy productions'
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
    # combine entropy productions as one structure
    ent = zeros(l,4)
    for i = 1:l
        ent[i,1] = prods[i,1] # A
        ent[i,2] = prods[i,3] # B
        ent[i,3] = ents[i,1] # A
        ent[i,4] = ents[i,3] # B
    end
    # combine entropy productions as one structure for Schlögl model
    Sent = zeros(l,4)
    for i = 1:l
        Sent[i,1] = Sprods[i,1] # A
        Sent[i,2] = Sprods[i,3] # B
        Sent[i,3] = Sents[i,1] # A
        Sent[i,4] = Sents[i,3] # B
    end
    # Line that sets up pyplot and basic sizings
    pyplot(dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12,tickfontsize=14)
    # Now a whole section for reading and plotting the stability data
    # first make structure to store data to
    acts = zeros(l,8)
    datayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2BD.csv"
            else
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2AD.csv"
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
                            acts[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
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
    # Now a whole section for reading and plotting the stability data for the Schlögl model
    # first make structure to store data to
    Sacts = zeros(l,8)
    Sdatayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])h2lD.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])l2hD.csv"
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
                            Sacts[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
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
    # plot change in action vs entropy produced
    xlab = L"ΔS_{A\rightarrow B} - ΔS_{B\rightarrow A}"
    ylab = L"\mathcal{A}_{B\rightarrow A} - \mathcal{A}_{A\rightarrow B}"
    plot(xlabel=xlab,ylabel=ylab,title=L"\Delta\mathcal{A}\;vs\;\Delta\Delta S")
    for i = 1:l
        if datayn[i] == true
            scatter!([acts[i,8]-acts[i,4]],[acts[i,2]-acts[i,6]],label="",color=1)
        end
    end
    for i = 1:l
        if Sdatayn[i] == true
            scatter!([Sacts[i,8]-Sacts[i,4]],[Sacts[i,2]-Sacts[i,6]],label="",color=2)
        end
    end
    # Now want to fit two lines, one for each model
    @. model(x,p) = p[1] + p[2]*x
    # Need to reduce data here to only include the relevant points
    xdataT = acts[datayn,8] .- acts[datayn,4]
    ydataT = acts[datayn,2] .- acts[datayn,6]
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xran = -18.0:1.0:52.0
    plot!(xran,model(xran,[yintT,slopT]),label="",color=1)
    # then same for Schlögl model
    xdataS = Sacts[Sdatayn,8] .- Sacts[Sdatayn,4]
    ydataS = Sacts[Sdatayn,2] .- Sacts[Sdatayn,6]
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    plot!(xran,model(xran,[yintS,slopS]),label="",color=2)
    savefig("../Results/DiffActvsDiffEntProd.png")
    # Now calculate Pearson correlation coefficient
    xbarT = sum(xdataT)/length(xdataT)
    ybarT = sum(ydataT)/length(ydataT)
    a = 0
    b = 0
    c = 0
    for i = 1:length(xdataT)
        a += (xdataT[i] - xbarT)*(ydataT[i] - ybarT)
        b += (xdataT[i] - xbarT)^2
        c += (ydataT[i] - ybarT)^2
    end
    r = a/sqrt(b*c)
    println("Correlation Toggle Switch Act vs Ent: $(r)")
    # repeat for Schlögl
    xbarS = sum(xdataS)/length(xdataS)
    ybarS = sum(ydataS)/length(ydataS)
    a = 0
    b = 0
    c = 0
    for i = 1:length(xdataS)
        a += (xdataS[i] - xbarS)*(ydataS[i] - ybarS)
        b += (xdataS[i] - xbarS)^2
        c += (ydataS[i] - ybarS)^2
    end
    r = a/sqrt(b*c)
    println("Correlation Schlögl Act vs Ent: $(r)")
    # Final section to read in probs and make a plot, this is necessarily rough at this point
    # Check there is a file of entropies to be read
    infile = "../Results/Fig3Data/$(ARGS[1])probs.csv"
    if ~isfile(infile)
        println("Error: No file of probabilities to be read.")
        return(nothing)
    end
    # now read in probabilities
    l = countlines(infile)
    w = 3
    probs = zeros(l,w)
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
                probs[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # now need to select inds based on whether probailities have been found yet
    ind = Array{Int64,1}(undef,0)
    for i = 1:l
        if probs[i,3] != 0.0
            ind = vcat(ind,i)
        end
    end
    # now read in probabilities for Schlögl model
    infile = "../Results/Fig3DataS/$(ARGS[2])probsS.csv"
    if ~isfile(infile)
        println("Error: No file of Schlögl probabilities to be read.")
        return(nothing)
    end
    l = countlines(infile)
    w = 3
    Sprobs = zeros(l,w)
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
                Sprobs[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # now need to select inds based on whether probabilities have been found yet
    Sind = Array{Int64,1}(undef,0)
    for i = 1:l
        if Sprobs[i,3] != 0.0
            Sind = vcat(Sind,i)
        end
    end
    plot(title=L"\ln{(\frac{p_{A}}{p_{B}})}\;vs\;\ln{(\frac{k_{B\rightarrow A}}{k_{A\rightarrow B}})}")
    for i = ind
        scatter!([probs[i,3]*(acts[i,2]-acts[i,6])],[log(probs[i,1]/probs[i,2])],label="",color=1)
    end
    for i = Sind
        scatter!([Sprobs[i,3]*(Sacts[i,2]-Sacts[i,6])],[log(Sprobs[i,1]/Sprobs[i,2])],label="",color=2)
    end
    x = -10.0:1.0:10.0
    plot!(x,x,label="",xlabel=L"\ln{(\frac{k_{B\rightarrow A}}{k_{A\rightarrow B}})}",ylabel=L"\ln{(\frac{p_{A}}{p_{B}})}")
    savefig("../Results/Linear.png")
    # now try to get log ratio of state probailities to plot against differences in entropy production
    ylab = L"\ln{(\frac{p_{A}}{p_{B}})}"
    xlab = L"\dot{S}_A - \dot{S}_B"
    plot(xlabel=xlab,ylabel=ylab,title="$(ylab) vs $(xlab)")
    for i = ind
        if datayn[i] == true
            scatter!([ent[i,3]-ent[i,4]],[log(probs[i,2]/probs[i,1])/probs[i,3]],label="",color=1)
        end
    end
    for i = Sind
        if Sdatayn[i] == true
            scatter!([Sent[i,3]-Sent[i,4]],[log(Sprobs[i,2]/Sprobs[i,1])/Sprobs[i,3]],label="",color=2)
        end
    end
    savefig("../Results/LogProbvsDiffEnt.png")
    # Now do graph of different
    plot(title=L"\Delta S_{A\rightarrow B} - \Delta S_{B\rightarrow A}",xlabel="Langevin EP (LDT)",ylabel="Exact EP (FT)")
    I = [2,3,11,13,19,32,48]
    for i = I
        infile1 = "../Results/Fig3Data/Traj/$(i)testA2BMast.csv"
        infile2 = "../Results/Fig3Data/Traj/$(i)testB2AMast.csv"
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
        # now read in productions and flows for foward and backwards paths
        infile1 = "../Results/Fig3Data/Traj/$(i)testA2Bpf.csv"
        infile2 = "../Results/Fig3Data/Traj/$(i)testB2Apf.csv"
        # Extract forward data
        w = 2
        pff = zeros(w)
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
                for i = 1:w
                    pff[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        pfb = zeros(w)
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
                for i = 1:w
                    pfb[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Now need to extract averages and standard deviations and plot against other data
        mf = mean(mastf)
        mb = mean(mastb)
        sdf = stdm(mastf,mf;corrected=true)/sqrt(ne-1)
        sdb = stdm(mastb,mb;corrected=true)/sqrt(ne-1)
        m = mf - mb
        sd = sqrt(sdf^2 + sdb^2)
        plot!([pff[1]-pfb[1]],[m],yerror=sd,label="")
    end
    # add linear lines
    x3 = -8.0:0.1:4.0
    plot!(x3,x3,label="")
    savefig("../Results/MultDiff.png")
    return(nothing)
end

@time main()
