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

# function to plot graphs of actions vs path entropies, path entropies vs prods, occupation probability from theory vs exact
# and occupation probability vs Schnakenberg entropy production
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
    xlab = L"ΔS^L_{B\rightarrow A} - ΔS^L_{A\rightarrow B}"
    ylab = L"\mathcal{A}_{A\rightarrow B} - \mathcal{A}_{B\rightarrow A}\;(1/\Omega)"
    plot(xlabel=xlab,ylabel=ylab,title="Path entropy production")
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
    # then same for Schlögl model
    xdataS = Sacts[Sdatayn,8] .- Sacts[Sdatayn,4]
    ydataS = Sacts[Sdatayn,2] .- Sacts[Sdatayn,6]
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    # Find confidence intervals of parameters and convert to usable form
    conT = confidence_interval(fitT, 0.05)
    vyintT = conT[1]
    vslopT = conT[2]
    conS = confidence_interval(fitS, 0.05)
    vyintS = conS[1]
    vslopS = conS[2]
    # Find deviating from best line and then plot them on raph as ribbon
    Tl = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[2],vslopT[2]]))
    Tu = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[1],vslopT[1]]))
    plot!(xran,model(xran,[yintT,slopT]),ribbon=(Tl,Tu),label="",color=1)
    # Same for Schlögl model
    Sl = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[2],vslopS[2]]))
    Su = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[1],vslopS[1]]))
    plot!(xran,model(xran,[yintS,slopS]),ribbon=(Sl,Su),label="",color=2)
    # Finally add points to entropy plot
    scatter!([xdataT],[ydataT],label="",color=1)
    scatter!([xdataS],[ydataS],label="",color=2)
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
    plot(title="Validity of FW action")
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
    # then same for Schlögl model
    xdataS = Sprobs[Sind,3].*(Sacts[Sind,2].-Sacts[Sind,6])
    ydataS = log.(Sprobs[Sind,1]./Sprobs[Sind,2])
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    # Find confidence intervals of parameters and convert to usable form
    conT = confidence_interval(fitT, 0.05)
    vyintT = conT[1]
    vslopT = conT[2]
    conS = confidence_interval(fitS, 0.05)
    vyintS = conS[1]
    vslopS = conS[2]
    # Find deviating from best line and then plot them on raph as ribbon
    Tl = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[2],vslopT[2]]))
    Tu = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[1],vslopT[1]]))
    plot!(xran,model(xran,[yintT,slopT]),ribbon=(Tl,Tu),label="",color=1)
    # Same for Schlögl model
    Sl = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[2],vslopS[2]]))
    Su = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[1],vslopS[1]]))
    plot!(xran,model(xran,[yintS,slopS]),ribbon=(Sl,Su),label="",color=2)
    # Finally add points to entropy plot
    scatter!([xdataT],[ydataT],label="",color=1)
    scatter!([xdataS],[ydataS],label="",color=2)
    plot!(xlabel=L"\mathcal{A}_{A\rightarrow B} - \mathcal{A}_{B\rightarrow A}\;(1/\Omega)")
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
    xlab = L"\dot{S}_A - \dot{S}_B\;(s^{-1})"
    plot(xlabel=xlab,ylabel=ylab,title="MaxEPP",ylims=(-45.0,22.5))
    xdataT = ent[ind,3].-ent[ind,4]
    ydataT = log.(probs[ind,1]./probs[ind,2])./probs[ind,3]
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xranT = -300.0:100.0:3800.0#1900.0
    # then same for Schlögl model
    xdataS = Sent[Sind,3].-Sent[Sind,4]
    ydataS = log.(Sprobs[Sind,1]./Sprobs[Sind,2])./Sprobs[Sind,3]
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    xranS = -300.0:100.0:4100.0
    # Find confidence intervals of parameters and convert to usable form
    conT = confidence_interval(fitT, 0.05)
    vyintT = conT[1]
    vslopT = conT[2]
    conS = confidence_interval(fitS, 0.05)
    vyintS = conS[1]
    vslopS = conS[2]
    # Find deviating from best line and then plot them on raph as ribbon
    Tl = abs.(model(xranT,[yintT,slopT]) .- model(xranT,[vyintT[2],vslopT[2]]))
    Tu = abs.(model(xranT,[yintT,slopT]) .- model(xranT,[vyintT[1],vslopT[1]]))
    plot!(xranT,model(xranT,[yintT,slopT]),ribbon=(Tl,Tu),label="",color=1)
    # Same for Schlögl model
    Sl = abs.(model(xranS,[yintS,slopS]) .- model(xranS,[vyintS[2],vslopS[2]]))
    Su = abs.(model(xranS,[yintS,slopS]) .- model(xranS,[vyintS[1],vslopS[1]]))
    plot!(xranS,model(xranS,[yintS,slopS]),ribbon=(Sl,Su),label="",color=2)
    # Now plot the data points on top
    for i = ind
        if datayn[i] == true
            scatter!([ent[i,3]-ent[i,4]],[log(probs[i,1]/probs[i,2])/probs[i,3]],label="",color=1)
        end
    end
    for i = Sind
        if Sdatayn[i] == true
            scatter!([Sent[i,3]-Sent[i,4]],[log(Sprobs[i,1]/Sprobs[i,2])/Sprobs[i,3]],label="",color=2)
        end
    end
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
    return(nothing)
end

# function to plot steady state entropy and diffusion matrix graphs
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
    # Check there is a file of magnitudes of diffusion matrices to be read
    infile = "../Results/Fig3Data/$(ARGS[1])D.csv"
    if ~isfile(infile)
        println("Error: No file of diffusion matrix magnitudes to be read for toggle.")
        return(nothing)
    end
    # now read in magnitudes of diffusion matrices
    l = countlines(infile)
    w = 3
    TDs = zeros(l,w)
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
                TDs[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Repeat for Schlögl model
    infile = "../Results/Fig3DataS/$(ARGS[2])DS.csv"
    if ~isfile(infile)
        println("Error: No file of diffusion matrix magnitudes to be read for Schlögl.")
        return(nothing)
    end
    # now read in magnitudes of diffusion matrices
    l = countlines(infile)
    w = 3
    SDs = zeros(l,w)
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
                SDs[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
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
    # Check there is a file of entropies to be read
    infile = "../Results/Fig3Data/$(ARGS[1])ent.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read for Toggle Switch.")
        return(nothing)
    end
    # now read in entropies
    l = countlines(infile)
    w = 2
    Tent = zeros(l,w)
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
                Tent[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Repeat for Schlögl model
    infile = "../Results/Fig3DataS/$(ARGS[2])entS.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read for Toggle Switch.")
        return(nothing)
    end
    # Now read in file of entropies
    l = countlines(infile)
    w = 2
    Sent = zeros(l,w)
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
                Sent[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Set basic pyplot settings
    pyplot(dpi=300,titlefontsize=20,guidefontsize=17,legendfontsize=15,tickfontsize=14)
    # Now try to make diffusion matrix plot
    scatter([acts[:,2].-acts[:,6]],[TDs[:,3]./TDs[:,1]],label="",color=1)
    scatter!([Sacts[:,2].-Sacts[:,6]],[SDs[:,3]./SDs[:,1]],label="",color=2)
    # Some Latex strings to insert
    xlab = L"\mathcal{A}_{A\rightarrow B} - \mathcal{A}_{B\rightarrow A}"
    ylab = L"mag(D_B)/mag(D_A)"
    titx = L"\Delta\mathcal{A}"
    plot!(xlabel=xlab,ylabel=ylab,title="$(ylab) vs $(titx)")
    savefig("../Results/DiffvsAct.png")
    # Some Latex strings to insert
    xlab = L"\mathcal{A}_{A\rightarrow B} - \mathcal{A}_{B\rightarrow A} (1/\Omega)"
    ylab = L"S_A - S_B"
    plot(xlabel=xlab,ylabel=ylab,title="Role of entropy")
    # Make linear model to fit to
    @. model(x, p) = p[1] + p[2]*x
    xdataT = acts[:,2].-acts[:,6]
    ydataT = Tent[:,1].-Tent[:,2]
    p0 = [0.0,1.0]
    fitT = curve_fit(model,xdataT,ydataT,p0)
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    xran = -20.0:1.0:28.0
    # then same for Schlögl model
    xdataS = Sacts[:,2].-Sacts[:,6]
    ydataS = Sent[:,1].-Sent[:,2]
    fitS = curve_fit(model,xdataS,ydataS,p0)
    yintS = coef(fitS)[1]
    slopS = coef(fitS)[2]
    # Find confidence intervals of parameters and convert to usable form
    conT = confidence_interval(fitT, 0.05)
    vyintT = conT[1]
    vslopT = conT[2]
    conS = confidence_interval(fitS, 0.05)
    vyintS = conS[1]
    vslopS = conS[2]
    # Find deviating from best line and then plot them on raph as ribbon
    Tl = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[2],vslopT[2]]))
    Tu = abs.(model(xran,[yintT,slopT]) .- model(xran,[vyintT[1],vslopT[1]]))
    plot!(xran,model(xran,[yintT,slopT]),ribbon=(Tl,Tu),label="",color=1)
    # Same for Schlögl model
    Sl = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[2],vslopS[2]]))
    Su = abs.(model(xran,[yintS,slopS]) .- model(xran,[vyintS[1],vslopS[1]]))
    plot!(xran,model(xran,[yintS,slopS]),ribbon=(Sl,Su),label="",color=2)
    # Finally add points to entropy plot
    scatter!([acts[:,2].-acts[:,6]],[Tent[:,1].-Tent[:,2]],label="",color=1)
    scatter!([Sacts[:,2].-Sacts[:,6]],[Sent[:,1].-Sent[:,2]],label="",color=2)
    savefig("../Results/EntvsAct.png")
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
    println("Correlation Toggle Entropy Difference vs Action Difference: $(r)")
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
    println("Correlation Schlögl Entropy Difference vs Action Difference: $(r)")
    return(nothing)
end

# A third function to investigate and plot the best and worst fitting parameter sets
# Only done for toggle switch so only need to provide one argument here
function third()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide a name for toggle switch inputs")
        return(nothing)
    end
    # First read in parameter sets
    infile = "../Results/Fig3Data/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 10
    ps = zeros(l,w)
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
    # Check file of steady states exists
    infile = "../Results/Fig3Data/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in steady states
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
    # Want to establish the worst and best fitting cases for toggle switch
    L = 10 # Number of points to find
    best = zeros(L,2)
    worst = zeros(L,2)
    δ = zeros(l)
    δ2 = zeros(l)
    r = zeros(l) # to store parameter ratios
    sr = zeros(l) # to store steady-state ratios
    fr = zeros(l) # to store f:r ratios
    for i = 1:l
        # Find divergence from expected relation
        δ[i] = abs(acts[i,8]-acts[i,4]-2*(acts[i,2]-acts[i,6]))
        sca = abs(acts[i,8]-acts[i,4]) # rescaling factor
        δ2[i] = δ[i]/sca
        sr[i] = (steads[i,3]+steads[i,4])/(steads[i,1]+steads[i,2]+steads[i,5]+steads[i,6])
        # Fill empty entries until vector full
        if i <= L
            if best[i,2] == 0.0 && worst[i,2] == 0.0
                best[i,1] = i
                worst[i,1] = i
                best[i,2] = δ2[i]
                worst[i,2] = δ2[i]
            end
        end
        # Slightly hacky reording when the two vectors are full
        if i == L+1
            bestp = sortperm(best[:,2])
            worstp = sortperm(worst[:,2],rev=true)
            worst = worst[worstp,:]
            best = best[bestp,:]
        end
        if i >= L+1
            # Nothing triggered yet
            tbest = false
            tworst = false
            for j = 1:L
                # Only run if not already triggered for this vector
                if best[j,2] > δ2[i] && tbest == false
                    for k = L:-1:j+1 # Use reverse for loop to rearrange elements
                        best[k,:] = best[k-1,:]
                    end
                    best[j,1] = i
                    best[j,2] = δ2[i]
                    tbest = true
                end
                if worst[j,2] < δ2[i] && tworst == false
                    for k = L:-1:j+1 # Use reverse for loop to rearrange elements
                        worst[k,:] = worst[k-1,:]
                    end
                    worst[j,1] = i
                    worst[j,2] = δ2[i]
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
    pyplot() # call pyplot
    # Now plot steady state ratios against rescaled divergences
    scatter(sr,δ2)
    savefig("../Results/BestWorst/SadvsDiv.png")
    return(nothing)
end

# little function to generate pyplot custom legends that I can then add to a bigger plot
function fourth()
    # call pyplot with basic settings
    pyplot(dpi=1000,titlefontsize=20,guidefontsize=17,legendfontsize=15,tickfontsize=14)
    x = 0.001:0.001:0.002
    plot(x,x,label=L"A{\rightarrow}B",color=:black,ticks=false,linestyle=:solid,ylim=(0.0,10000.0))
    plot!(x,x,label=L"B{\rightarrow}A",color=:black,linestyle=:dash,border=:none,xlim=(0.0,10000.0))
    savefig("../Results/CustomLegend.png")
    return(nothing)
end

@time first()
@time second()
