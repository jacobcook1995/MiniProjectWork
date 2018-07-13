#!/usr/bin/env julia
# gill4.jl
# A script to efficently perform a gillespie simulation of the full 4 species model
# This script generates  histograms of the probability distribution from gillespie simulation

using Roots
using Plots
import GR

# function to find the zeros of the function
function nullcline(F::Float64,r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,
                    q::Float64,kmin::Float64,qmin::Float64,Ne::Float64)
    g(x) = (K + kmin)*(q/k)*((r + f*((F - K*x)/Q)^2)/(r + f*x^2))*x - (qmin + Q)*(F - K*x)/Q
    three = false
    n = 0
    As = []
    while three == false
        As = fzeros(g, 0, 2*F/K, order = 1)
        n = length(As)
        if n == 3
            three = true
        end
    end
    Bs = zeros(n)
    Ss = zeros(n)
    Ws = zeros(n)
    for i = 1:n
        Bs[i] = (F - K*As[i])/Q
        Ss[i] = (1/(k*r))*(r + f*((F - K*As[i])/Q)^2)*(K + kmin)*As[i]
        Ws[i] = Ne - As[i] - Bs[i] - Ss[i]
    end
    ss1 = [ As[1]; Bs[1]; Ss[1]; Ws[1] ]
    sad = [ As[2]; Bs[2]; Ss[2]; Ws[2] ]
    ss2 = [ As[3]; Bs[3]; Ss[3]; Ws[3] ]
    print("$(ss1)\n")
    print("$(sad)\n")
    print("$(ss2)\n")
    return (ss1,sad,ss2)
end

# function to construct the rates
function rates(A::Int64,B::Int64,S::Int64,W::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,F::Float64)
    rates = [r*k*S/(r+f*B*(B-1)), kmin*A, K*A, Kmin*W, r*q*S/(r+f*A*(A-1)), qmin*B, Q*B, Qmin*W, F]
    return(rates)
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates:::Array{Float64,1}vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars[1] += 1 # A produced, S consumed
        vars[3] -= 1
    elseif r < rs[1] + rs[2]
        vars[1] -= 1
        vars[3] += 1 # A degraded producing an S
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays to a W
        vars[4] += 1
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1
        vars[4] -= 1 # W reforms as an A
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5]
        vars[2] += 1
        vars[3] -= 1 # B produced, S consumed
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6]
        vars[2] -= 1
        vars[3] += 1 # B degraded producing an S
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7]
        vars[2] -= 1
        vars[4] += 1 # B decays to a W
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7] + rs[8]
        vars[2] += 1
        vars[4] -= 1 # W reforms to a B
    else
        vars[3] += 1
        vars[4] -= 1 # W removed and S supplied
    end
    return(vars)
end

# function to actually run a gillepsie simulation
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,F::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1},
                    maxA::Int64,maxB::Int64,maxS::Int64,maxW::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,4,2)
    vars[:,2] = star
    # need to think about the histogram sizes here
    histA = zeros(3*maxA)
    histB = zeros(3*maxB)
    histS = zeros(3*maxS)
    histW = zeros(3*maxW)
    for i = 1:noits
        vars[:,1] = vars[:,2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[1,1],vars[2,1],vars[3,1],vars[4,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,F)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[:,2] = step(rs,vars[:,1])
        # add to histogram
        histA[vars[1,1]+1] += times[2] - times[1]
        histB[vars[2,1]+1] += times[2] - times[1]
        histS[vars[3,1]+1] += times[2] - times[1]
        histW[vars[4,1]+1] += times[2] - times[1]
    end
    histA = histA/times[2]
    histB = histB/times[2]
    histS = histS/times[2]
    histW = histW/times[2]
    println("Gillespie Done!")
    return(histA,histB,histS,histW)
end

# main function
function main()
    # General parameters
    Ω = 2 # Ω = 1, gMAP parameterisation
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/(Ω^2) # Promoter switching
    r = 10.0
    F = 10.0*Ω
    Kmin = 10.0^-20 # remains neligable though
    Qmin = 10.0^-20
    Ne = 150.0*Ω # number of elements in the system

    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(F,r,f,K,Q,k,q,kmin,qmin,Ne)
    # round star so that it becomes a vector of integers
    star2 = fill(0,4)
    mid2 = fill(0,4)
    fin2 = fill(0,4)
    for i = 1:4
        star2[i] = round(Int64,star1[i])
        mid2[i] = round(Int64,mid1[i])
        fin2[i] = round(Int64,fin1[i])
    end
    # find maximum value of each variable
    maxA = maximum([star2[1], mid2[1], fin2[1]])
    maxB = maximum([star2[2], mid2[2], fin2[2]])
    maxS = maximum([star2[3], mid2[3], fin2[3]])
    maxW = maximum([star2[4], mid2[4], fin2[4]])
    # now run gillespie
    noits = 500000000
    # maybe randomise the starting point somewhat
    # also maybe remove the vars
    histA, histB, histS, histW = gillespie(K,k,Q,q,kmin,qmin,f,F,r,Kmin,Qmin,noits,mid2,maxA,maxB,maxS,maxW)
    lisA = collect(0:(3*maxA-1))
    lisB = collect(0:(3*maxB-1))
    lisS = collect(0:(3*maxS-1))
    lisW = collect(0:(3*maxW-1))
    bar(lisA,histA)
    savefig("../Results/MeGraph.png")
    bar(lisB,histB)
    savefig("../Results/MeGraph2.png")
    bar(lisS,histS)
    savefig("../Results/MeGraph3.png")
    bar(lisW,histW)
    savefig("../Results/MeGraph4.png")
end

@time main()
