#!/usr/bin/env julia
# langsch.jl
# A script that impliments an overdamped langevin simulation of the two models
# this is done to generate a probability distribution which will then be used
# with the method of Tomé (2006) to find dissapted power
#
# Author: Jacob Cook
# Date: November 2018
using Roots
using DifferentialEquations
using Plots
using StatsBase

# ps = [ k1, K1, k2, K2, B, V, D ]
# total force
function f!(ps::Array{Float64,1},f::Float64,pos::Float64)
    f = ps[1] - ps[2]*pos + ps[4]*ps[5]*(pos^2) - ps[3]*(pos^3)
    return(f)
end
# Conservative force
function fc!(ps::Array{Float64,1},f::Float64,pos::Float64)
    f = ps[1] - ps[2]*pos + ps[4]*ps[5]*(pos^2) - ps[3]*(pos^3)
    return(f)
end
# Dissaptive force
function fd!(ps::Array{Float64,1},f::Float64,pos::Float64)
    f = 0.0
    return(f)
end
# noise
function fii!(ps::Array{Float64,1},f::Float64,pos::Float64)
    f = -ps[2] + 2*ps[4]*ps[5]*pos - 3*ps[3]*(pos^2)
    return(f)
end
# D
function D!(ps::Array{Float64,1},D::Float64,pos::Float64)
    D = ps[7]
    return(D)
end
function g(du::Array{Float64,1},u::Array{Float64,1},ps::Array{Float64,1},t::Float64)
    du[1] = ps[1] - ps[2]*u[1] + ps[4]*ps[5]*(u[1]^2) - ps[3]*(u[1]^3)
    return(du)
end

# noise
function σ(du::Array{Float64,},u::Array{Float64,1},ps::Array{Float64,1},t::Float64)
    D = ps[7]
    V = ps[6]
    noi = sqrt(2.0*D)/V
    du[1] = noi
    return(du)
end

# function to generate steady states
function states(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    # write out equation to be solved
    f(x) = k1 - K1*x - k2*(x^3) + K2*B*(x^2)
    three = false
    Xs = [ 0.0, 0.0, 0.0 ]
    while three == false
        X = fzeros(f, 0, 10*V, no_pts = 1000)
        if length(X) == 3
            three = true
            Xs = X
        end
    end
    ss1 = [Xs[1]]
    sad = [Xs[2]]
    ss2 = [Xs[3]]
    return(ss1,sad,ss2)
end

# function to provide normalised parameters for the gillespie simulation
function paras(V::Float64)
    # just gonna hard code some in for the moment
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    B = 4.0
    D = 1.0
    # refine constants by V
    k1 = k1*V
    k2 = k2/(V^2)
    K2 = K2/(V)
    return(k1,K1,k2,K2,B,D)
end

function generate()
    filename = "../Results/LangData/dataS.csv"
    V = 20.0
    k1, K1, k2, K2, B, D = paras(V)
    pars = [ k1, K1, k2, K2, B, V, D ]
    if isfile(filename)
        line = open(readlines, filename)[end]
        comma = 0
        L = length(line)
        i = 1
        for j = 1:L
            if line[j] == ','
                comma = j
                i += 1
            end
        end
        X = parse(Float64, line[1:(comma - 1)])
        ss1 = [X]
    else
        _, _, ss1 = states(k1,K1,k2,K2,B,V)
        # ss1 = V*ss1
    end
    println(ss1)
    # make a new function for specific volume and parameters
    prob_sde_lorenz = SDEProblem(g,σ,ss1,(0.0,10000.0),pars)
    sol = solve(prob_sde_lorenz)
    # step to convert from awkward Array{Array} form
    Xs = hcat(sol.u...)'
    data = zeros(length(Xs),2)
    # make data array to output
    for i = 1:length(Xs)
        data[i,1] = Xs[i]
        if i == 1
            data[i,2] = (sol.t[i+1] - sol.t[i])/2
        elseif i == length(Xs)
            data[i,2] = (sol.t[i] - sol.t[i-1])/2
        else
            data[i,2] = (sol.t[i+1] - sol.t[i-1])/2
        end
    end
    # now write function to print out this data
    if isfile(filename)
        # if file already exists should add to bottom of it
        out_file = open(filename, "a")
        for i = 1:size(data,1)
            line = "$(data[i,1]),$(data[i,2])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        out_file = open(filename, "w")
        for i = 1:size(data,1)
            line = "$(data[i,1]),$(data[i,2])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    return(nothing)
end

function analysis()
    filename = "../Results/LangData/dataS.csv"
    V = 20.0
    k1, K1, k2, K2, B, D = paras(V)
    pars = [ k1, K1, k2, K2, B, V, D ]
    # set up fast way to read file
    no_lins = countlines(filename)
    data = Array{Float64,2}(undef,no_lins,2)
    lincount = 1
    open(filename, "r") do in_file
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            i = 1
            for j = 1:L
                if line[j] == ','
                    comma = j
                    i += 1
                end
            end
            X = parse(Float64, line[1:(comma - 1)])
            dt = parse(Float64, line[(comma + 1):L])
            data[lincount,:] = [ X dt ]
            lincount += 1
        end
    end
    data[:,2] = data[:,2]/sum(data[:,2])
    # make histogram weighted by times
    points = data[:,1]
    hist = fit(Histogram,points,weights(data[:,2]),closed=:left,nbins=500)
    Xrange = hist.edges[1]
    lenX = Xrange.len
    # create arrays to store values of interest
    ΦC = zeros(lenX-1)
    ΦD = zeros(lenX-1)
    f = 0.0
    fc = 0.0
    fd = 0.0
    fii = 0.0
    Ds = 0.0
    pos = 0.0 # vector to provide posistions
    # now run for loop to finds f's etc
    for i = 1:lenX-1
        # want midpoint here
        pos = (Xrange[i] + Xrange[i+1])/2
        f = f!(pars,f,pos)
        fc = fc!(pars,fc,pos)
        fd = fd!(pars,fd,pos)
        fii = fii!(pars,fii,pos)
        Ds = D!(pars,Ds,pos)
        ΦC[i] = (fc[1]*f[1]/Ds[1] + sum(fii))*hist.weights[i]
        ΦD[i] = (fd[1]*f[1]/Ds[1])*hist.weights[i]
    end
    histC = deepcopy(hist)
    histD = deepcopy(hist)
    histP = deepcopy(hist)
    histC.weights = V*ΦC
    histD.weights = V*ΦD
    histP.weights = V*(ΦC.+ΦD)
    plot(hist)
    savefig("../Results/test.png")
    plot(histC)
    savefig("../Results/test1.png")
    plot(histD)
    savefig("../Results/test2.png")
    plot(histP)
    savefig("../Results/test3.png")
    println(sum(V*ΦC))
    println(sum(V*ΦD))
    return(nothing)
end

# @time generate()
for i = 1:10
    @time generate()
end
@time analysis()
