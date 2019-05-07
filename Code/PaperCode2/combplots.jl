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
    pyplot(dpi=300,titlefontsize=20,guidefontsize=17,legendfontsize=15,tickfontsize=14)
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
    xlab = L"ΔS^L_{A\rightarrow B} - ΔS^L_{B\rightarrow A}"
    ylab = L"\mathcal{A}_{B\rightarrow A} - \mathcal{A}_{A\rightarrow B}\;(1/\Omega)"
    plot(xlabel=xlab,ylabel=ylab,title=L"\Delta\mathcal{A}\;vs\;\Delta\Delta S^L")#,bgcolor=:transparent,fgcolor=:black)
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
    # Want to establish the worst and best fitting cases for toggle switch
    best = zeros(10,2)
    worst = zeros(10,2)
    for i = 1:l
        # Find divergence from expected relation
        δ = abs(acts[i,8]-acts[i,4]-2*(acts[i,2]-acts[i,6]))
        # δ = abs(acts[i,8]-acts[i,4]-2*(acts[i,2]-acts[i,6]))/(abs(acts[i,8]-acts[i,4])+abs(2*(acts[i,2]-acts[i,6])))
        # Fill empty entries until vector full
        if i <= 10
            if best[i,2] == 0.0 && worst[i,2] == 0.0
                best[i,1] = i
                worst[i,1] = i
                best[i,2] = δ
                worst[i,2] = δ
            end
        end
        # Slightly hacky reording when the two vectors are full
        if i == 11
            bestp = sortperm(best[:,2])
            worstp = sortperm(worst[:,2],rev=true)
            worst = worst[worstp,:]
            best = best[bestp,:]
        end
        if i >= 11
            # Nothing triggered yet
            tbest = false
            tworst = false
            for j = 1:10
                # Only run if not already triggered for this vector
                if best[j,2] > δ && tbest == false
                    for k = 10:-1:j+1 # Use reverse for loop to rearrange elements
                        best[k,:] = best[k-1,:]
                    end
                    best[j,1] = i
                    best[j,2] = δ
                    tbest = true
                end
                if worst[j,2] < δ && tworst == false
                    for k = 10:-1:j+1 # Use reverse for loop to rearrange elements
                        worst[k,:] = worst[k-1,:]
                    end
                    worst[j,1] = i
                    worst[j,2] = δ
                    tworst = true
                end
                # break if both are triggered
                if tbest == true && tworst == true
                    break
                end
            end
        end
    end
    # Then output reduced best and worst indices
    outfileb = "../Results/Fig3Data/Best.csv"
    outfilew = "../Results/Fig3Data/Worst.csv"
    # open files for writing
    out_fileb = open(outfileb, "w")
    for j = 1:size(best,1)
        line = "$(best[j,1]),$(best[j,2])\n"
        write(out_fileb, line)
    end
    close(out_fileb)
    out_filew = open(outfilew, "w")
    for j = 1:size(worst,1)
        line = "$(worst[j,1]),$(worst[j,2])\n"
        write(out_filew, line)
    end
    close(out_filew)
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
    plot!(xran,model(xran,[yintS,slopS]),label="",color=2,top_margin=8mm)
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
    plot(title=L"\ln{\left(\frac{P_{A}}{P_{B}}\right)}\;vs\;\ln{\left(\frac{k_{B\rightarrow A}}{k_{A\rightarrow B}}\right)}")#,bgcolor=:transparent,fgcolor=:black)
    for i = ind
        scatter!([probs[i,3]*(acts[i,2]-acts[i,6])],[log(probs[i,1]/probs[i,2])],label="",color=1)
    end
    for i = Sind
        scatter!([Sprobs[i,3]*(Sacts[i,2]-Sacts[i,6])],[log(Sprobs[i,1]/Sprobs[i,2])],label="",color=2)
    end
    xdataT = probs[ind,3].*(acts[ind,2].-acts[ind,6])
    ydataT = log.(probs[ind,1]./probs[ind,2])
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xran = -10.0:1.0:13.0
    plot!(xran,model(xran,[yintT,slopT]),label="",color=1)
    # then same for Schlögl model
    xdataS = Sprobs[Sind,3].*(Sacts[Sind,2].-Sacts[Sind,6])
    ydataS = log.(Sprobs[Sind,1]./Sprobs[Sind,2])
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    plot!(xran,model(xran,[yintS,slopS]),label="",color=2,xlabel=L"\ln{\left(\frac{k_{B\rightarrow A}}{k_{A\rightarrow B}}\right)}")
    plot!(xran,xran,label="",ylabel=L"\ln{\left(\frac{P_{A}}{P_{B}}\right)}",color=:black,style=:dash)
    savefig("../Results/Linear.png")
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
    println("Correlation Toggle Switch Theory vs Simulation: $(r)")
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
    println("Correlation Schlögl Theory vs Simulation: $(r)")
    # now try to get log ratio of state probailities to plot against differences in entropy production
    ylab = L"\ln{\left(\frac{P_{A}}{P_{B}}\right)}"
    xlab1 = L"\dot{S}_A - \dot{S}_B"
    xlab2 = L"\dot{S}_A - \dot{S}_B\;(s^{-1})"
    plot(xlabel=xlab2,ylabel=ylab,title="$(ylab) vs $(xlab1)")#,bgcolor=:transparent,fgcolor=:black)
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
    xdataT = ent[ind,3].-ent[ind,4]
    ydataT = log.(probs[ind,2]./probs[ind,1])./probs[ind,3]
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xran = -300.0:100.0:1900.0
    plot!(xran,model(xran,[yintT,slopT]),label="",color=1)
    # then same for Schlögl model
    xdataS = Sent[Sind,3].-Sent[Sind,4]
    ydataS = log.(Sprobs[Sind,2]./Sprobs[Sind,1])./Sprobs[Sind,3]
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    xran = -300.0:100.0:4100.0
    plot!(xran,model(xran,[yintS,slopS]),label="",color=2)
    savefig("../Results/LogProbvsDiffEnt.png")
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
    println("Correlation Toggle Switch Occ vs Ent Prod: $(r)")
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
    println("Correlation Schlögl Occ vs Ent Prod: $(r)")
    # Read in steady states and use to calculate distances
    infile = "../Results/Fig3Data/$(ARGS[1])stead.csv"
    l = countlines(infile)
    w = 6
    steads = zeros(l,w)
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
    # now make δsad
    δsad = zeros(l,2)
    for i = 1:l
        δsad[i,1] = sqrt((steads[i,1]-steads[i,3])^2 + (steads[i,2]-steads[i,4])^2)
        δsad[i,2] = sqrt((steads[i,5]-steads[i,3])^2 + (steads[i,6]-steads[i,4])^2)
    end
    # Same for Schlögl states
    infile = "../Results/Fig3DataS/$(ARGS[2])steadS.csv"
    l = countlines(infile)
    w = 3
    steadsS = zeros(l,w)
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
                steadsS[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # now make δsad
    δsadS = zeros(l,2)
    for i = 1:l
        δsadS[i,1] = abs(steadsS[i,1]-steadsS[i,2])
        δsadS[i,2] = abs(steadsS[i,3]-steadsS[i,2])
    end
    # Now do graph of different
    xlab = L"ΔS"
    ylab = L"\dot{S}"
    mag1 = 10^-2
    Lmag1 = L"10^2"
    p1 = plot(title=L"\Delta S_{A\rightarrow B} - \Delta S_{B\rightarrow A}",guidefontsize=15)#,bgcolor=:transparent,fgcolor=:black)
    plot!(p1,xlabel="Langevin EP (LDT) ($Lmag1)",ylabel="Exact EP (FT) ($Lmag1)")
    # p2 = plot(xlabel=xlab,ylabel=ylab,title=L"Δ\;\dot{S}\;vs\;\Delta\Delta S")
    # p3 = plot(xlabel=xlab,ylabel=ylab,title=L"Δ\;\dot{S}\;vs\;\Delta\Delta S")
    I = collect(1:100)
    K = 0 # start counter
    # make arrays to store data
    pos = zeros(length(I))
    ms = zeros(length(I))
    sds = zeros(length(I))
    # pos2 = zeros(length(I))
    # save all x and y values
    for i = I
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
        # now read in productions and flows for foward and backwards paths
        infile1 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2Bpf.csv"
        infile2 = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2Apf.csv"
        # Extract forward data
        w = 4
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
        plot!(p1,mag1*[pff[1]-pfb[1]],mag1*[m],yerror=mag1*sd,label="",markercolor=1,markerstrokecolor=:black)
        # scatter!(p2,[pff[3]-pfb[3]],[ents[i,1] - ents[i,3]],label="",color=1)
        # scatter!(p3,[pff[3]-pfb[3]],[δsad[i,1]*ents[i,1] - δsad[i,2]*ents[i,3]],label="",color=1)
        # Save data for use outside
        K += 1 # increment index
        pos[K] = pff[1]-pfb[1]
        ms[K] = m
        sds[K] = sd
        # pos2[K] = pff[3]-pfb[3]
    end
    # Abandoning theory here and just finding the best fit
    xdataT = pos
    ydataT = ms
    weigT = (sds.^-2)
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,weigT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xran = -55.0:5.0:195.0
    plot!(p1,mag1*xran,mag1*model(xran,[yintT,slopT]),label="",color=1)
    # xdataT = pos2
    # ydataT = ents[:,1] .- ents[:,3]
    # fitT = curve_fit(model,xdataT,ydataT,p0)
    # yintT = coef(fitT)[1]
    # slopT = coef(fitT)[2]
    # xran = -45.0:5.0:10.0
    # plot!(p2,xran,model(xran,[yintT,slopT]),label="",color=1)
    # xdataT = pos2
    # ydataT = δsad[:,1].*ents[:,1] .- δsad[:,2].*ents[:,3]
    # fitT = curve_fit(model,xdataT,ydataT,p0)
    # yintT = coef(fitT)[1]
    # slopT = coef(fitT)[2]
    # xran = -45.0:5.0:10.0
    # plot!(p3,xran,model(xran,[yintT,slopT]),label="",color=1)
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
    println("Correlation Toggle Switch EntProd vs TrajEntProd: $(r)")
    # And could do a weighted correlation
    wxbarT = sum(weigT.*xdataT)/(length(xdataT)*sum(weigT))
    wybarT = sum(weigT.*ydataT)/(length(ydataT)*sum(weigT))
    wcovxy = sum(weigT.*(xdataT.-wxbarT).*(ydataT.-wybarT))/sum(weigT)
    wcovxx = sum(weigT.*(xdataT.-wxbarT).*(xdataT.-wxbarT))/sum(weigT)
    wcovyy = sum(weigT.*(ydataT.-wybarT).*(ydataT.-wybarT))/sum(weigT)
    r = wcovxy/sqrt(wcovxx*wcovyy)
    println("Weighted Correlation Toggle Switch Mast vs Langevin: $(r)")
    # The same is now done for the Schlögl model
    # make arrays to store data
    pos = zeros(length(I))
    ms = zeros(length(I))
    sds = zeros(length(I))
    # pos2 = zeros(length(I))
    K = 0 # restart counter
    # save all x and y values
    for i = I
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
        # now read in productions and flows for foward and backwards paths
        infile1 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])h2lpf.csv"
        infile2 = "../Results/Fig3DataS/Traj/$(i)$(ARGS[2])l2hpf.csv"
        # Extract forward data
        w = 4
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
        plot!(p1,mag1*[pff[1]-pfb[1]],mag1*[m],yerror=mag1*sd,label="",markercolor=2,markerstrokecolor=:black)
        # scatter!(p2,[pff[3]-pfb[3]],[Sents[i,1] - Sents[i,3]],label="",color=2)
        # scatter!(p3,[pff[3]-pfb[3]],[δsadS[i,1]*Sents[i,1] - δsadS[i,2]*Sents[i,3]],label="",color=2)
        # Save data for use outside
        K += 1 # increment index
        pos[K] = pff[1]-pfb[1]
        ms[K] = m
        sds[K] = sd
        # pos2[K] = pff[3]-pfb[3]
    end
    # Abandoning theory here and just finding the best fit
    xdataS = pos
    ydataS = ms
    weigS = (sds.^-2)
    p0 = [0.0,1.0]
    fitS = curve_fit(model,xdataS,ydataS,weigS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    xran = -55.0:5.0:195.0
    plot!(p1,mag1*xran,mag1*model(xran,[yintS,slopS]),label="",color=2)
    # xdataS = pos2
    # ydataS = Sents[:,1] .- Sents[:,3]
    # p0 = [0.0,1.0]
    # fitS = curve_fit(model,xdataS,ydataS,p0)
    # yintS = coef(fitS)[1]
    # slopS = coef(fitS)[2]
    # xran = -45.0:5.0:10.0
    # plot!(p2,xran,model(xran,[yintS,slopS]),label="",color=2)
    # xdataS = pos2
    # ydataS = δsadS[:,1].*Sents[:,1] .- δsadS[:,2].*Sents[:,3]
    # p0 = [0.0,1.0]
    # fitS = curve_fit(model,xdataS,ydataS,p0)
    # yintS = coef(fitS)[1]
    # slopS = coef(fitS)[2]
    # xran = -45.0:5.0:10.0
    # plot!(p3,xran,model(xran,[yintS,slopS]),label="",color=2)
    savefig(p1,"../Results/MultDiff.png")
    # savefig(p2,"../Results/SwitchvsProd.png")
    # savefig(p3,"../Results/SwitchvsProdtest.png")
    # Now calculate Pearson correlation coefficient
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
    println("Correlation Schlögl Entprod vs TrajEntProd: $(r)")
    # And could do a weighted correlation
    wxbarS = sum(weigS.*xdataS)/(length(xdataS)*sum(weigS))
    wybarS = sum(weigS.*ydataS)/(length(ydataS)*sum(weigS))
    wcovxy = sum(weigS.*(xdataS.-wxbarS).*(ydataS.-wybarS))/sum(weigS)
    wcovxx = sum(weigS.*(xdataS.-wxbarS).*(xdataS.-wxbarS))/sum(weigS)
    wcovyy = sum(weigS.*(ydataS.-wybarS).*(ydataS.-wybarS))/sum(weigS)
    r = wcovxy/sqrt(wcovxx*wcovyy)
    println("Weighted Correlation Schlögl Mast vs Langevin: $(r)")
    return(nothing)
end

@time main()
