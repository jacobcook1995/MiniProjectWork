#!/usr/bin/env julia
# gilltau.jl
# A script to efficently perform a gillespie simulation of the reduced 2 species model
# This script generates  histograms of the probability distribution from gillespie simulation

using Roots
using Plots
import GR

# finction to find start end and saddle points
function nullcline(r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64)
    a = 2
    b = 2
    A1(x) = k*r/(K*(r+f*x^a))
    A2(x) = (r/f*(q/(Q*x)-1))^(1/b)
    g(x) = k*r/(K*(r+f*x^a)) - (r/f*(q/(Q*x)-1))^(1/b) #A1(x) - A2(x)
    xs = fzeros(g, 0, q/Q)
    ss1 = [A1(xs[1]); xs[1]]
    sad = [A1(xs[2]); xs[2]]
    ss2 = [A1(xs[3]); xs[3]]
    print(ss1)
    print("\n")
    print(sad)
    print("\n")
    print(ss2)
    print("\n")
    return (ss1,sad,ss2)
end

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,diffs::Bool=false)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    if diffs == false
        return(rates)
    else
        dA = [ 1, -1, -1, 1, 0, 0, 0, 0]
        dB = [ 0, 0, 0, 0, 1, -1, -1, 1]
        return(rates,dA,dB)
    end
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
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,f::Float64,
                    r::Float64,Kmin::Float64,Qmin::Float64,tf::Float64,star::Array{Int64,1},Ω::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,2,2)
    vars[:,2] = star
    τ1 = 0
    fin = false
    while fin == false
        vars[:,1] = vars[:,2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[:,2] = step(rs,vars[:,1])
        if vars[2,2] < star[2]
            τ1 += τ
        elseif vars[2,2] == star[2]
            τ1 += 0.5*τ
        end
        if times[2] >= tf
            fin = true
        end
    end
    τ1 = τ1/times[2]
    return(τ1)
end

# main function
function main()
    # General parameters
    Ω = 150
    K = 10.0
    k = K*Ω # steady state for A=k/K=1
    Q = 1.0
    q = Q*Ω
    kmin = 10.0^-20 # set all too 10.0^-20 for now
    Kmin = 10.0^-20
    qmin = 10.0^-20
    Qmin = 10.0^-20
    f = 1000.0/(Ω^2) # Promoter switching
    r = 10.0

    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(r,f,K,Q,k,q,kmin,qmin)
    # round star so that it becomes a vector of integers
    star2 = fill(0,2)
    mid2 = fill(0,2)
    fin2 = fill(0,2)
    for i = 1:2
        star2[i] = round(Int64,star1[i])
        mid2[i] = round(Int64,mid1[i])
        fin2[i] = round(Int64,fin1[i])
    end
    # now run gillespie
    noruns = 50000
    tf = 0.02*(1)
    τ = zeros(noruns)
    for i = 1:noruns
        τ[i] = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,tf,fin2,Ω)
    end
    # construct my own histogram
    lis = collect(0:0.001:1)
    lis2 = collect(0.005:0.001:0.995)
    l2 = length(lis2)
    τhist = zeros(l2)
    st = round(Int64,l2/2)
    for i = 1:noruns
        j = st
        found = false
        while found == false
            if τ[i] >= lis[j]
                if j == l2 || τ[i] < lis[j+1]
                    τhist[j] += 1
                    found = true
                else
                    j += 1
                end
            else
                j -= 1
            end
        end
    end
    bar(lis2,τhist)
    savefig("../Results/Hist$(tf)B.png")
    return(nothing)
end

@time main()
