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
pygui(:qt5)
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
    # change this so that it uses less memory
    ABf = zeros(1001,2)
    ABb = zeros(1001,2)
    forw = false
    back = false
    ABf0 = copy(star)
    ABf[1,:] = ABf0
    ABb0 = copy(fin)
    ABb[1,:] = ABb0
    # for i = 1:noits
    #     # calculate rates
    #     rs = rates(AB[1],AB[2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
    #     # do gillepsie step
    #     AB = step(rs,AB)
    #     # add to vectors
    #     ABs[i+1,:] = AB
    # end
    println("Gillespie Done!")
    return(ABf,ABb)
end

# main function
function main()
    # General parameters
    Ω = 150
    k = 26.097758831774257
    kmin = 2.179419011537317
    q = 25.772181639137496
    qmin = 1.4148644973536413
    K = 8.24009582574649
    Kmin = 0.6004846952009286
    Q = 6.96376978178525
    Qmin = 0.1657402116759403
    r = 3875.552664673383
    f = 9285.774344574023
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
    # switch to PyPlot
    pyplot()
    # DEfine relevant latex strings
    LatS = latexstring("\\Omega = $(Ω)")
    # Then plot graphs
    plot(ABf[:,2]/Ω,ABf[:,1]/Ω,label="to high A",title="2D Toggle Switch",legendtitle=LatS)
    plot!(ABb[:,2]/Ω,ABb[:,1]/Ω,label="to high B",xlabel="A",ylabel="B")
    scatter!([star2[2]/Ω], [star2[1]/Ω], seriescolor = :green, label="")
    scatter!([mid2[2]/Ω], [mid2[1]/Ω], seriescolor = :orange, label="")
    scatter!([fin2[2]/Ω], [fin2[1]/Ω], seriescolor = :red, label="")
    savefig("../Results/switch.png")
    return(nothing)
end

@time main()
