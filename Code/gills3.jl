#!/usr/bin/env julia
# gills3.jl
# A script to efficently perform a gillespie simulation of the full 3 species schlogl model
# This script generates  histograms of the probability distribution from gillespie simulation

using Roots
using Plots
import GR

# finction to find start end and saddle points
function nullcline(k1::Float64,K1::Float64,k2::Float64,K2::Float64,F::Float64,N::Float64)
    # write out equation to be solved
    three = false
    Xs = [ ]
    n = 0
    f(x) = -k2*x^3 + 4.0*K2*x^2 - K1*x + 1.0*k1
    while three == false
        Xs = fzeros(f, 0, N, order = 1)
        println(Xs)
        n = length(Xs)
        if n == 2#3
            three = true
        end
    end
    As = zeros(n)
    Bs = zeros(n)
    for i = 1:n
        As[i] = (K1*Xs[i] + F)/k1
        Bs[i] = (k2*Xs[i]^3 - F)/(K2*Xs[i]^2)
    end
    ss1 = [ Xs[1], As[1], Bs[1] ]
    println(ss1)
    sad = [ Xs[2], As[2], Bs[2] ]
    println(sad)
    ss2 = [ Xs[3], As[3], Bs[3] ]
    println(ss2)
    return(ss1,sad,ss2)
end

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,V::Int64)
    rates = [ k1*V, K1*X, k2*X*(X-1)*(X-2)/(V*V), K2*B*X*(X-1)/V]
    return(rates)
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars += 1 # X produced from A (not simulated)
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels to A (not simulated)
    elseif r < rs[1] + rs[2] + rs[3]
        vars -= 1 # X decays to B (not simulated)
    else
        vars += 1 # X regenerated from B (not simulated)
    end
    return(vars)
end

# function to actually run a gillepsie simulation
function gillespie(k1::Float64,K1::Float64,k2::Float64,K2::Float64,noits::Int64,star::Int64,Ω::Int64,maxX::Int64,maxA::Int64,maxB::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,3,2)
    vars[:,2] = star
    histX = zeros(4*maxX)
    histA = zeros(4*maxA)
    histB = zeros(4*maxB)
    for i = 1:noits
        vars[:,1] = vars[:,2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[:,1],k1,K1,k2,K2,Ω)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[:,2] = step(rs,vars[:,1])
        # add to histogram
        histX[vars[1,1]+1] += times[2] - times[1]
        histA[vars[2,1]+1] += times[2] - times[1]
        histB[vars[3,1]+1] += times[2] - times[1]
    end
    histX = histX/times[2]
    histA = histA/times[2]
    histB = histB/times[2]
    println("Gillespie Done!")
    return(histX,histA,histB)
end

# main function
function main()
    # General parameters
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    Ω = 20
    F = -1.0
    N = 100.0


    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(k1,K1,k2,K2,F,N)
    # round star so that it becomes a vector of integers
    star2 = fill(0,3)
    mid2 = fill(0,3)
    fin2 = fill(0,3)
    for i = 1:3
        star2[i] = round(Int64,star1[i])
        mid2[i] = round(Int64,mid1[i])
        fin2[i] = round(Int64,fin1[i])
    end

    # now run gillespie
    noits = 500000000
    maxX = maximum([star2[1], mid2[1], fin2[1]])
    maxA = maximum([star2[2], mid2[2], fin2[2]])
    maxB = maximum([star2[3], mid2[3], fin2[3]])
    histX, histA, histB = gillespie(k1,K1,k2,K2,noits,mid2,Ω,maxX,maxA,maxB)
    lisX = collect(0:(4*maxX-1))
    lisA = collect(0:(4*maxA-1))
    lisB = collect(0:(4*maxB-1))
    bar(lisX,histX)
    savefig("../Results/MeGraph.png")
    bar(lisA,histA)
    savefig("../Results/MeGraph2.png")
    bar(lisB,histB)
    savefig("../Results/MeGraph3.png")
    return(nothing)
end

@time main()
