#!/usr/bin/env julia
# lang2spec.jl
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

# ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω ]
# total force
function f!(ps::Array{Float64,1},f::Array{Float64,1},pos::Array{Float64,1})
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*pos[2]^2) + ps[6] - (ps[2]+ps[5])*pos[1]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*pos[1]^2) + ps[8] - (ps[4]+ps[7])*pos[2]
    return(f)
end
# Conservative force
function fc!(ps::Array{Float64,1},f::Array{Float64,1},pos::Array{Float64,1})
    f[1] = -(ps[2]+ps[5])*pos[1]
    f[2] = -(ps[4]+ps[7])*pos[2]
    return(f)
end
# Dissaptive force
function fd!(ps::Array{Float64,1},f::Array{Float64,1},pos::Array{Float64,1})
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*pos[2]^2) + ps[6]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*pos[1]^2) + ps[8]
    return(f)
end
# noise
function fii!(ps::Array{Float64,1},f::Array{Float64,1},pos::Array{Float64,1})
    f[1] = -(ps[2]+ps[5])
    f[2] = -(ps[4]+ps[7])
    return(f)
end
# D
function D!(ps::Array{Float64,1},D::Array{Float64,1},pos::Array{Float64,1})
    D[1] = ps[12]
    D[2] = ps[12]
    return(D)
end
function g(du::Array{Float64,1},u::Array{Float64,1},ps::Array{Float64,1},t::Float64)
    du[1] = ps[1]*ps[9]/(ps[9]+ps[10]*u[2]^2) + ps[6] - (ps[2]+ps[5])*u[1]
    du[2] = ps[3]*ps[9]/(ps[9]+ps[10]*u[1]^2) + ps[8] - (ps[4]+ps[7])*u[2]
    return(du)
end

# noise
function σ(du::Array{Float64,1},u::Array{Float64,1},ps::Array{Float64,1},t::Float64)
    D = ps[12]
    Ω = ps[11]
    noi = sqrt(2.0*D)/Ω
    du[1] = noi
    du[2] = noi
    return(du)
end

# function to generate steady states
function states(k::Float64,kmin::Float64,q::Float64,qmin::Float64,K::Float64,Kmin::Float64,
                Q::Float64,Qmin::Float64,r::Float64,f::Float64,Ω::Int64)
    # unnormalise the parameters as this is at Langevin level
    k = k/Ω
    Kmin = Kmin/Ω
    q = q/Ω
    Qmin = Qmin/Ω
    f = f*(Ω^2)
    # define two nullcline equations
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
    return(ss1,sad,ss2)
end

# function to provide normalised parameters for the gillespie simulation
function paras()
    # just gonna hard code some in for the moment
    # EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT
    Ω = 200
    D = 1.0
    k = 10.0
    kmin = 0.1
    q = 5.0
    qmin = 0.1
    K = 1.0
    Kmin = 0.01
    Q = 0.5
    Qmin = 0.01
    r = 1000.0
    f = 10000.0
    # then normalise appropriatly
    k = k*Ω
    Kmin = Kmin*Ω
    q = q*Ω
    Qmin = Qmin*Ω
    f = f/(Ω^2)
    return(k,kmin,q,qmin,K,Kmin,Q,Qmin,r,f,Ω,D)
end

function generate()
    filename = "../Results/LangData/data.csv"
    k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω, D = paras(Ω)
    pars = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω, D ]
    if isfile(filename)
        line = open(readlines, filename)[end]
        comma = fill(0,2)
        L = length(line)
        i = 1
        for j = 1:L
            if line[j] == ','
                comma[i] = j
                i += 1
            end
        end
        A = parse(Float64, line[1:(comma[1] - 1)])
        B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
        ss1 = [ A, B ]
    else
        ss1, _, _ = states(k,kmin,q,qmin,K,Kmin,Q,Qmin,r,f,Ω)
        ss1 = Ω*ss1
    end
    println(ss1)
    # make a new function for specific volume and parameters
    prob_sde_lorenz = SDEProblem(g,σ,ss1,(0.0,1000000.0),pars)
    sol = solve(prob_sde_lorenz)
    # step to convert from awkward Array{Array} form
    AB = hcat(sol.u...)'
    data = zeros(size(AB,1),3)
    # make data array to output
    for i = 1:size(AB,1)
        data[i,1] = AB[i,1]
        data[i,2] = AB[i,2]
        if i == 1
            data[i,3] = (sol.t[i+1] - sol.t[i])/2
        elseif i == size(AB,1)
            data[i,3] = (sol.t[i] - sol.t[i-1])/2
        else
            data[i,3] = (sol.t[i+1] - sol.t[i-1])/2
        end
    end
    # now write function to print out this data
    if isfile(filename)
        # if file already exists should add to bottom of it
        out_file = open(filename, "a")
        for i = 1:size(data,1)
            line = "$(data[i,1]),$(data[i,2]),$(data[i,3])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        out_file = open(filename, "w")
        for i = 1:size(data,1)
            line = "$(data[i,1]),$(data[i,2]),$(data[i,3])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    return(nothing)
end

function analysis()
    filename = "../Results/LangData/data.csv"
    k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω, D = paras()
    pars = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω, D ]
    # set up fast way to read file
    no_lins = countlines(filename)
    data = Array{Float64,2}(undef,no_lins,3)
    lincount = 1
    open(filename, "r") do in_file
        for line in eachline(in_file)
            # parse line by finding commas
            comma = fill(0,2)
            L = length(line)
            i = 1
            for j = 1:L
                if line[j] == ','
                    comma[i] = j
                    i += 1
                end
            end
            A = parse(Float64, line[1:(comma[1] - 1)])
            B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
            dt = parse(Float64, line[(comma[2] + 1):L])
            data[lincount,:] = [ A B dt ]
            lincount += 1
        end
    end
    data[:,3] = data[:,3]/sum(data[:,3])
    # make histogram weighted by times
    points = tuple(data[:,1],data[:,2])
    println("READ IN DATA")
    hist = fit(Histogram,points,weights(data[:,3]),closed=:left,nbins=500)
    println("MADE HISTOGRAM")
    Arange = hist.edges[1]
    Brange = hist.edges[2]
    lenA = Arange.len
    lenB = Brange.len
    # create arrays to store values of interest
    ΦC = zeros(lenA-1,lenB-1)
    ΦD = zeros(lenA-1,lenB-1)
    ΦS = zeros(lenA-1,lenB-1)
    ΦS2 = zeros(lenA-1,lenB-1)
    f = zeros(2)
    fc = zeros(2)
    fd = zeros(2)
    fii = zeros(2)
    Ds = zeros(2)
    pos = zeros(2) # vector to provide posistions
    # now run for loop to finds f's etc
    for j = 1:lenB-1
        pos[2] = (Brange[j] + Brange[j+1])/2
        for i = 1:lenA-1
            # want midpoint here
            pos[1] = (Arange[i] + Arange[i+1])/2
            f = f!(pars,f,pos)
            fc = fc!(pars,fc,pos)
            fd = fd!(pars,fd,pos)
            fii = fii!(pars,fii,pos)
            Ds = D!(pars,Ds,pos)
            ΦC[i,j] = (fc[1]*f[1]/Ds[1] + fc[2]*f[2]/Ds[2] + sum(fii))*hist.weights[i,j]
            ΦD[i,j] = (fd[1]*f[1]/Ds[1] + fd[2]*f[2]/Ds[2])*hist.weights[i,j]
            ΦS[i,j] = (f[1]*f[1]/Ds[1] + f[2]*f[2]/Ds[2])*hist.weights[i,j]
        end
    end
    histC = deepcopy(hist)
    histD = deepcopy(hist)
    histP = deepcopy(hist)
    histS = deepcopy(hist)
    histC.weights = Ω*ΦC
    histD.weights = Ω*ΦD
    histP.weights = Ω*(ΦC.+ΦD)
    histS.weights = Ω*(ΦS)
    plot(hist)
    savefig("../Results/test.png")
    plot(histC)
    savefig("../Results/test1.png")
    plot(histD)
    savefig("../Results/test2.png")
    plot(histP)
    savefig("../Results/test3.png")
    plot(histS)
    savefig("../Results/test4.png")
    println(sum(Ω*ΦC))
    println(sum(Ω*ΦD))
    return(nothing)
end
# for i = 1:5
#     @time generate()
# end
@time analysis()
