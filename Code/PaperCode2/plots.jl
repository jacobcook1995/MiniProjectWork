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
using StatsPlots
using DataFrames
using Statistics
using Plots.PlotMeasures

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
    # Check there is a file of entropies to be read
    infile = "../Results/Fig3Data/$(ARGS[1])ent.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read.")
        return(nothing)
    end
    # now read in entropies
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
    # Check there is a file of magnitudes of diffusion matrices to be read
    infile = "../Results/Fig3Data/$(ARGS[1])D.csv"
    if ~isfile(infile)
        println("Error: No file of entropies to be read.")
        return(nothing)
    end
    # now read in magnitudes of diffusion matrices
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
    pyplot(dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12,tickfontsize=14)
    # scatter([ent[:,3]],[Ds[:,1]],label="")
    # scatter!([ent[:,4]],[Ds[:,3]],label="")
    # plot!(xlabel=L"\dot{S}",ylabel=L"mag(D)",title="Magnitude of D vs Entropy Production")
    # savefig("../Results/DvsEntProd.png")
    # scatter([ent[:,3].-ent[:,4]],[Ds[:,1]./Ds[:,3]],label="")
    # plot!(xlabel=L"\dot{S}_{h}-\dot{S}_{l}",ylabel=L"mag(D_{h})/mag(D_{l})",title="Ratio of D vs Steady State Ent Prod Diff")
    # savefig("../Results/DiffDvsEntProd.png")
    # scatter([shent[:,1].-shent[:,2]],[Ds[:,1]./Ds[:,3]],label="")
    # plot!(xlabel=L"S_{h}-S_{l}",ylabel=L"mag(D_{h})/mag(D_{l})",title="Ratio of D vs Steady State Ent Diff")
    # savefig("../Results/DiffDvsEnt.png")
    # scatter([ent[:,1].-ent[:,2]],[ent[:,3].-ent[:,4]],label="")
    # plot!(xlabel=L"\dot{S}_{h,prod}-\dot{S}_{l,prod}",ylabel=L"\dot{S}_h - \dot{S}_l",title="Diff in Entropy vs Production")
    # savefig("../Results/DifferProds.png")
    # # Plot entropies against Schnakenberg entropy production
    # scatter([ent[:,3]],[shent[:,1]],label="")
    # scatter!([ent[:,4]],[shent[:,2]],label="")
    # plot!(xlabel=L"\dot{S}",ylabel=L"S",title="Ent vs Ent Prod")
    # savefig("../Results/EntvsSchnak.png")
    scatter([ent[:,3].-ent[:,4]],[shent[:,1].-shent[:,2]],label="",ylim=(-10,10))
    plot!(xlabel=L"\dot{S}_h-\dot{S}_l",ylabel=L"S_h - S_l",title=L"\Delta S\;vs\;\Delta\dot{S}",top_margin=8.0mm)
    savefig("../Results/DiffEnt.png")
    # # plot scatter graph of the wrong points
    # plot(xlabel="Ent Prod Rate Term",ylabel="Ent Prod Rate Traj",title="Mismatched points")
    # for i = 1:length(wrong)
    #     scatter!([ent[wrong[i],1],ent[wrong[i],2]], [ent[wrong[i],5],ent[wrong[i],6]],label="")
    # end
    # savefig("../Results/WrongTrajvsTerms.png")
    # # then plot scatter graph
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
    # plot change in action vs entropy produced
    xlab = L"ΔS_{h\rightarrow l} - ΔS_{l\rightarrow h}"
    ylab = L"\mathcal{A}_{h\rightarrow l} - \mathcal{A}_{l\rightarrow h}"
    plot(xlabel=xlab,ylabel=ylab,title=L"\Delta\mathcal{A}\;vs\;\Delta\Delta S")
    for i = 1:l
        if datayn[i] == true
            scatter!([acts[i,4]-acts[i,8]],[acts[i,2]-acts[i,6]],label="",color=1)
        end
    end
    savefig("../Results/DiffActvsDiffEntProd.png")
    # plot(xlabel="Entropy Produced",ylabel="Action",title="Action vs Entropy Produced")
    # for i = 1:l
    #     if datayn[i] == true
    #         scatter!([acts[i,4],acts[i,8]],[acts[i,2],acts[i,6]],label="")
    #     end
    # end
    # savefig("../Results/ActvsEntProd.png")
    # plot(xlabel="Approx Action",ylabel="Action",title="Action vs Approx Action")
    # for i = 1:l
    #     if datayn[i] == true
    #         scatter!([acts[i,3],acts[i,7]],[acts[i,2],acts[i,6]],label="")
    #     end
    # end
    # savefig("../Results/ActvsActapprox.png")
    # # now try to get log ratio of rate of switching to plot against differences in entropy production
    # lab = L"\ln{\frac{k_{l\rightarrow h}}{k_{h\rightarrow l}}}"
    # labx = L"\dot{S}_h - \dot{S}_l"
    # plot(xlabel=labx,ylabel=lab,title="$(lab) vs $(labx)")
    # for i = 1:l
    #     if datayn[i] == true
    #         scatter!([ent[i,3]-ent[i,4]],[acts[i,6]-acts[i,2]],label="",color=1)
    #     end
    # end
    # savefig("../Results/LogStabvsDiffEnt.png")
    # now plot differnce in entropy production along path vs differnce at steady states
    # lab = L"ΔS_{h\rightarrow l} - ΔS_{l\rightarrow h}"
    # plot(xlabel=L"\dot{S}_h - \dot{S}_l",ylabel=lab,title="Diff in Path Ent Prod vs Steady State")
    # for i = 1:l
    #     if datayn[i] == true
    #         scatter!([ent[i,3]-ent[i,4]],[acts[i,4]-acts[i,8]],label="")
    #     end
    # end
    # savefig("../Results/PathvsState.png")
    # Final section to read in probs and make a plot, this is necessarily rough at this point
    # Check there is a file of entropies to be read
    infile = "../Results/Fig3Data/$(ARGS[1])probs.csv"
    if ~isfile(infile)
        println("Error: No file of probabilities to be read.")
        return(nothing)
    end
    # now read in entropies
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
    p1 = plot(title=L"\ln{\frac{p_{h}}{p_{l}}}\;vs\;\ln{\frac{k_{l\rightarrow h}}{k_{h\rightarrow l}}}")
    # p2 = plot(title="Residuals")
    for i = ind
        scatter!(p1,[probs[i,3]*(acts[i,2]-acts[i,6])],[log(probs[i,1]/probs[i,2])],label="",color=1)
    #     scatter!(p2,[probs[i,3]*(acts[i,2]-acts[i,6])-log(probs[i,1]/probs[i,2])],[0.0],label="")
    end
    x = -7.5:0.1:4.5
    plot!(p1,x,x,label="",xlabel=L"\ln{\frac{k_{l\rightarrow h}}{k_{h\rightarrow l}}}",ylabel=L"\ln{\frac{p_{h}}{p_{l}}}")
    savefig(p1,"../Results/Linear.png")
    # savefig(p2,"../Results/resid.png")
    # now try to get log ratio of state probailities to plot against differences in entropy production
    ylab = L"\ln{\frac{p_{h}}{p_{l}}}"
    xlab = L"\dot{S}_h - \dot{S}_l"
    plot(xlabel=xlab,ylabel=ylab,title="$(ylab) vs $(xlab)")
    for i = ind
        if datayn[i] == true
            scatter!([ent[i,3]-ent[i,4]],[log(probs[i,2]/probs[i,1])/probs[i,3]],label="",color=1)
        end
    end
    savefig("../Results/LogProbvsDiffEnt.png")
    # # now try to get log ratio of rate of switching to plot against differences in entropy production
    # lab = L"\ln{\frac{k_{l\rightarrow h}}{k_{h\rightarrow l}}}"
    # plot(xlabel=L"\dot{S}_h - \dot{S}_l",ylabel=lab,title="Log ratio of switching vs Ent Prod Diff")
    # for i = ind
    #     if datayn[i] == true
    #         scatter!([ent[i,3]-ent[i,4]],[acts[i,6]-acts[i,2]],label="")
    #     end
    # end
    # savefig("../Results/LogStabvsDiffEnt.png")
    # # Now read in my two written data csvs for trajctories
    # infile1 = "../Results/Fig3Data/Traj/2testA2BMas.csv"
    # infile2 = "../Results/Fig3Data/Traj//2testB2AMas.csv"
    # w = 101
    # mastf = zeros(100,10)
    # open(infile1, "r") do in_file
    #     # Use a for loop to process the rows in the input file one-by-one
    #     k = 1
    #     for line in eachline(in_file)
    #         # parse line by finding commas
    #         L = length(line)
    #         comma = fill(0,w+1)
    #         j = 1
    #         for i = 1:L
    #             if line[i] == ','
    #                 j += 1
    #                 comma[j] = i
    #             end
    #         end
    #         comma[end] = L+1
    #         for i = 1:w-1
    #             mastf[i,k] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
    #         end
    #         k += 1
    #     end
    # end
    # mastb = zeros(100,10)
    # open(infile2, "r") do in_file
    #     # Use a for loop to process the rows in the input file one-by-one
    #     k = 1
    #     for line in eachline(in_file)
    #         # parse line by finding commas
    #         L = length(line)
    #         comma = fill(0,w+1)
    #         j = 1
    #         for i = 1:L
    #             if line[i] == ','
    #                 j += 1
    #                 comma[j] = i
    #             end
    #         end
    #         comma[end] = L+1
    #         for i = 1:w-1
    #             mastb[i,k] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
    #         end
    #         k += 1
    #     end
    # end
    # # do box plot of forward data first
    # df = DataFrame(a=mastf[:,1],b=mastf[:,2],c=mastf[:,3],d=mastf[:,4],e=mastf[:,5],f=mastf[:,6],g=mastf[:,7],h=mastf[:,8],i=mastf[:,9],j=mastf[:,10])
    # @df df boxplot([:a :b :c :d :e :f :g :h :i :j],label=["" "" "" "" "" "" "" "" "" ""])
    # hline!([acts[2,4]],label = L"\Delta S",xlabel=L"\Omega",ylabel="Entropy Produced",title="Ent Prod A→B")
    # savefig("../Results/A2B.png")
    # df = DataFrame(a=mastb[:,1],b=mastb[:,2],c=mastb[:,3],d=mastb[:,4],e=mastb[:,5],f=mastb[:,6],g=mastb[:,7],h=mastb[:,8],i=mastb[:,9],j=mastb[:,10])
    # @df df boxplot([:a :b :c :d :e :f :g :h :i :j],label=["" "" "" "" "" "" "" "" "" ""])
    # hline!([acts[2,8]],label = L"\Delta S",xlabel=L"\Omega",ylabel="Entropy Produced",title="Ent Prod B→A")
    # savefig("../Results/B2A.png")
    # a = (sum(mastf[:,1])-sum(mastb[:,1]))/100
    # b = (sum(mastf[:,2])-sum(mastb[:,2]))/100
    # c = (sum(mastf[:,3])-sum(mastb[:,3]))/100
    # d = (sum(mastf[:,4])-sum(mastb[:,4]))/100
    # e = (sum(mastf[:,5])-sum(mastb[:,5]))/100
    # f = (sum(mastf[:,6])-sum(mastb[:,6]))/100
    # g = (sum(mastf[:,7])-sum(mastb[:,7]))/100
    # h = (sum(mastf[:,8])-sum(mastb[:,8]))/100
    # i = (sum(mastf[:,9])-sum(mastb[:,9]))/100
    # j = (sum(mastf[:,10])-sum(mastb[:,10]))/100
    # Ωss = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
    # means = [a b c d e f g h i j]
    # errors = zeros(10)
    # for i = 1:10
    #     mf = mean(mastf[:,i])
    #     mb = mean(mastb[:,i])
    #     sdf = stdm(mastf[:,i],mf;corrected=true)/sqrt(w-1)
    #     sdb = stdm(mastb,mb;corrected=true)/sqrt(w-1)
    #     errors[i] = sqrt(sdf^2 + sdb^2)
    # end
    # plot(Ωss,means,yerror=errors,label="",xlabel=L"\Omega",ylabel=L"\Delta S_{A\rightarrow B} - \Delta S_{B\rightarrow A}")
    # hline!([acts[2,4]-acts[2,8]],label="",title="Diff Ent Prod")
    # savefig("../Results/Diff.png")
    # Now do graph of different
    p1 = plot(title="A→B",xlabel="Entropy Prod Langevin",ylabel="Entropy Prod Master")
    p2 = plot(title="B→A",xlabel="Entropy Prod Langevin",ylabel="Entropy Prod Master")
    p3 = plot(title=L"\Delta S_{A\rightarrow B} - \Delta S_{B\rightarrow A}",xlabel="E_P",ylabel="Master")
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
        plot!(p1,[pff[1]],[mf],yerror=sdf,label="")
        plot!(p2,[pfb[1]],[mb],yerror=sdb,label="")
        plot!(p3,[pff[1]-pfb[1]],[m],yerror=sd,label="")
    end
    # add linear lines
    x1 = -0.0:1.0:800.0
    plot!(p1,x1,x1,label="")
    x2 = -0.0:1.0:800.0
    plot!(p2,x2,x2,label="")
    x3 = -8.0:0.1:4.0
    plot!(p3,x3,x3,label="")
    savefig(p1,"../Results/MultA2B.png")
    savefig(p2,"../Results/MultB2A.png")
    savefig(p3,"../Results/MultDiff.png")
    return(nothing)
end

@time main()
