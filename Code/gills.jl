#!/usr/bin/env julia
# gills.jl
# A script to efficently perform a gillespie simulation of the reduced 1 species schlogl model
# This script generates  histograms of the probability distribution from gillespie simulation

using Roots
using Plots
import GR

# finction to find start end and saddle points
function nullcline(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    # write out equation to be solved
    f(x) = k1 - K1*x - k2*(x^3) + K2*B*(x^2)
    three = false
    Xs = [ 0.0, 0.0, 0.0 ]
    while three == false
        X = fzeros(f, 0, 10, no_pts = 1000)
        if length(X) == 3
            three = true
            Xs = X
        end
    end
    Xs = Xs*V
    ss1 = Xs[1]
    sad = Xs[2]
    ss2 = Xs[3]
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
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
function step(rates::Array{Float64,1},vars::Int64)
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
function gillespie(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,noits::Int64,star::Int64,Ω::Int64,maxX::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,2)
    vars[2] = star
    histX = zeros(4*maxX)
    for i = 1:noits
        vars[1] = vars[2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[1],k1,K1,k2,K2,B,Ω)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[2] = step(rs,vars[1])
        # add to histogram
        histX[vars[1]+1] += times[2] - times[1]
    end
    histX = histX/times[2]
    println("Gillespie Done!")
    return(histX)
end

# main function
function main()
    # General parameters
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    B = 4.0
    Ω = 20

    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(k1,K1,k2,K2,B,Ω)
    # round star so that it becomes a vector of integers
    star2 = round(Int64,star1)
    mid2 = round(Int64,mid1)
    fin2 = round(Int64,fin1)

    # now run gillespie
    noits = 500000000
    maxX = maximum([star2, mid2, fin2])
    histX = gillespie(k1,K1,k2,K2,B,noits,mid2,Ω,maxX)
    lis = collect(0:(4*maxX-1))
    bar(lis,histX)
    savefig("../Results/MeGraph.png")
    return(nothing)
end

@time main()
