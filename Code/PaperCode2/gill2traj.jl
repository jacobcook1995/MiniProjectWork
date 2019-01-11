#!/usr/bin/env julia
# gill2traj.jl
# A script to make nice plots of a two species gillespie simulation showing a forward
# and backwards trajectory between the two states
#
# Author: Jacob Cook
# Date: January 2019

using Roots
using Plots
using LaTeXStrings
using PyCall
import PyPlot

# function to find start end and saddle points
function nullcline(r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    a = 2
    b = 2
    A1(x) = real(sqrt(complex((r/f)*(q/((qmin+Q)*x - Qmin) - 1))))
    A2(x) = (1/(kmin+K))*((k*r)/(r + f*x^2) + Kmin)
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    # then setup loop to find three solutions to these equations
    while three == false
        bs1 = fzeros(g, 0.0, 0.1) # catches zeros about the origin
        bs2 = fzeros(g, 0.1, 2*q/Q) # the upper bound here is potentially problematic
        bs = vcat(bs1,bs2)
        n = length(bs)
        gs = 0
        bad = zeros(Int64,0)
        for i = 1:n
            # check if this is actual solution and not artifact of the numerical method
            gs = g(bs[i])
            tol = 1.0e-14
            if gs >= 0 + tol || gs <= 0 - tol
                bad = append!(bad,i)
            end
        end
        if length(bad) != 0
            n = n - length(bad)
            bs = deleteat!(bs, bad)
        end
        if n == 3
            three = true
        end
    end
    # now make the stationary points
    ss1 = [ A1(bs[1]), bs[1] ]
    sad = [ A1(bs[2]), bs[2] ]
    ss2 = [ A1(bs[3]), bs[3] ]
    println("The three stationary points are:")
    println(ss1)
    println(sad)
    println(ss2)
    flush(stdout) # needed to ensure output prints when using nohup
    return(ss1,sad,ss2)
end

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars[1] += 1 # A produced
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1 # A regenerated
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5]
        vars[2] += 1 # B produced
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6]
        vars[2] -= 1 # B unravels
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7]
        vars[2] -= 1 # B decays
    else
        vars[2] += 1 # B regenerated
    end
    return(vars)
end

# function to actually run a gillepsie simulation
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,star::Array{Int64,1},
                    fin::Array{Int64,1},Ω::Int64)
    # rescale constants appropraiatly
    Kmin = Kmin*Ω
    Qmin = Qmin*Ω
    k = k*Ω
    q = q*Ω
    f = f/(Ω^2)
    # set up arrarys to store paths
    L = 100000000
    ABf = zeros(2,L)
    ABb = zeros(2,L)
    # first do forward transistion
    forw = false
    ABf0 = copy(star)
    ABf[:,1] = ABf0
    j = forj = 0
    t1 = time_ns()
    while forw == false
        j += 1
        # calculate rates
        rs = rates(ABf0[1],ABf0[2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # do gillepsie step
        ABf0 = step(rs,ABf0)
        # add to vectors
        ABf[:,j] = ABf0
        # stopping condition
        if ABf0[1] <= fin[1] && ABf0[2] >= fin[2]
            forw = true
            forj = j
        elseif j >= L # overflow condition
            j = 1
            ABf[:,1] = ABf[:,end]
        end
    end
    t2 = time_ns()
    println("Forward Gillespie Done!")
    println("Time Elapsed: $((t2-t1)*10.0^-9)")
    flush(stdout)
    # and then backwards transisition
    back = false
    ABb0 = copy(fin)
    ABb[:,1] = ABb0
    j = bacj = 0
    t1 = time_ns()
    while back == false
        j += 1
        # calculate rates
        rs = rates(ABb0[1],ABb0[2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # do gillepsie step
        ABb0 = step(rs,ABb0)
        # add to vectors
        ABb[:,j] = ABb0
        # stopping condition
        if ABb0[1] >= star[1] && ABb0[2] <= star[2]
            back = true
            bacj = j
        elseif j >= L # overflow condition
            j = 1
            ABb[:,1] = ABb[:,end]
        end
    end
    t2 = time_ns()
    println("Backward Gillespie Done!")
    println("Time Elapsed: $((t2-t1)*10.0^-9)")
    flush(stdout)
    # now remove invalid data from array
    ABf2 = ABf[:,1:forj]
    ABb2 = ABb[:,1:bacj]
    fori = baci = 1
    # now search for sensible starting point for trajectories
    for i = forj:(-1):2
        if ABf2[1,i] >= star[1] && ABf2[2,i] <= star[2]
            fori = i
            break
        end
    end
    for i = bacj:(-1):2
        if ABb2[1,i] <= fin[1] && ABb2[2,i] >= fin[2]
            baci = i
            break
        end
    end
    # remove uneccesary earlier trajectory
    ABf3 = ABf2[:,fori:end]
    ABb3 = ABb2[:,baci:end]
    return(ABf3,ABb3)
end

# main function
function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # General parameters
    Ω = 80
    k = 20.315460339280257
    kmin = 0.19687208653639934
    q = 12.294415636568056
    qmin = 0.37021152930020407
    K = 7.792600383813772
    Kmin = 0.7093969626035892
    Q = 4.669640667899318
    Qmin = 0.3211077205202387
    r = 1516.1520996175293
    f = 7513.713695501023
    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    # round star so that it becomes a vector of integers
    star2 = fill(0,2)
    mid2 = fill(0,2)
    fin2 = fill(0,2)
    for i = 1:2
        star2[i] = round(Int64,star1[i]*Ω)
        mid2[i] = round(Int64,mid1[i]*Ω)
        fin2[i] = round(Int64,fin1[i]*Ω)
    end
    # now run gillespie
    ABf, ABb = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,star2,fin2,Ω)
    # First write out parameters
    outfile = "../Results/Fig1Data/$(ARGS[1])paras.csv"
    out_file = open(outfile, "w")
    ps = [k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f , Ω, star2[1], star2[2], mid2[1], mid2[2], fin2[1], fin2[2]]
    line = ""
    for j = 1:length(ps)
        line *= "$(ps[j]),"
    end
    # remove surplus ,
    line = line[1:end-1]
    # add a new line
    line *= "\n"
    write(out_file,line)
    close(outfile)
    # Then write out forward trajectory
    outfile1 = "../Results/Fig1Data/$(ARGS[1])for.csv"
    out_file1 = open(outfile1, "w")
    for j = 1:size(ABf,2)
        line = "$(ABf[1,j]),$(ABf[2,j])\n"
        write(out_file1, line)
    end
    close(outfile1)
    # Then write out backward trajectory
    outfile2 = "../Results/Fig1Data/$(ARGS[1])bac.csv"
    out_file2 = open(outfile2, "w")
    for j = 1:size(ABb,2)
        line = "$(ABb[1,j]),$(ABb[2,j])\n"
        write(out_file2, line)
    end
    close(outfile2)
    return(nothing)
end

# Seperate plotting function
function plots()
    println("Compiled, Starting plotting function.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to use to find input.")
        return(nothing)
    end
    infile = "../Results/Fig1Data/$(ARGS[1])paras.csv"
    w = 17
    ps = zeros(w)
    open(infile, "r") do in_file
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
                ps[i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
        end
    end
    Ω = ps[11]
    star2 = [ps[12], ps[13]]
    mid2 = [ps[14], ps[15]]
    fin2 = [ps[16], ps[17]]
    # Now read in forward path
    infile = "../Results/Fig1Data/$(ARGS[1])for.csv"
    w = 2
    len = countlines(infile)
    ABf = zeros(len,w)
    open(infile, "r") do in_file
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
                ABf[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now read in backward path
    infile = "../Results/Fig1Data/$(ARGS[1])bac.csv"
    w = 2
    len = countlines(infile)
    ABb = zeros(len,w)
    open(infile, "r") do in_file
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
                ABb[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # switch to PyPlot
    pyplot()
    L1 = L"\bullet\rightarrow\circ"
    L2 = L"\circ\rightarrow\bullet"
    # Then plot graphs
    plot(ABf[:,2],ABf[:,1],label=L1,title="2D Toggle Switch",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=15)
    plot!(ABb[:,2],ABb[:,1],label=L2,xlabel="Copy number A",ylabel="Copy number B")
    scatter!([star2[2]],[star2[1]],markersize=6,markercolor=:black,label="")
    scatter!([mid2[2]],[mid2[1]],markersize=5,markercolor=:black,markershape=:x,label="")
    scatter!([fin2[2]],[fin2[1]],markersize=6,markercolor=:white,label="")
    savefig("../Results/switch.png")
end

# @time main()
@time plots()
