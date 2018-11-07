#!/usr/bin/env julia
# trajentgill.jl
# A script to find the entropy production of the gillespie simulation of the schlögl model
#
# Author: Jacob Cook
# Date: November 2018

using Roots
using Plots
import GR

# Now construct the three relevant vectors of equations
function f!(F::Float64,X::Float64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    F = k1*V - K1*X - k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return F
end

# Then construct the necessary matrices
function D!(D::Float64,X::Float64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    D = k1*V + K1*X + k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return D
end

# finction to find start end and saddle points
function nullcline(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,high2low::Bool,V::Int64)
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
    sad = Xs[2]
    if high2low == true
        ss1 = Xs[3]
        ss2 = Xs[1]
    else
        ss1 = Xs[1]
        ss2 = Xs[3]
    end
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

# Schnakenberg formula for entropy Production
function Sch(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,pos::Float64)
    # fluxes
    r1 = k1*V
    rmin1 = K1*pos
    r2 = k2*(pos^3)/(V^2)
    rmin2 = K2*B*(pos^2)/V
    # net fluxes
    F1 = r1 - rmin1
    F2 = r2 - rmin2
    # affinities
    A1 = log(r1/rmin1)
    A2 = log(r2/rmin2)
    # entropy productions
    S1 = F1*A1
    S2 = F2*A2
    S = S1 + S2
    return(S)
end

# Schnakenberg formula for reduced form of the entropy Production
function Schr(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,pos::Float64)
    # fluxes
    r1 = k1*V
    r2 = K1*pos
    r3 = k2*(pos^3)/(V^2)
    r4 = K2*B*(pos^2)/V
    # net flux
    F = r1 + r4 - r2 - r3
    # affinity
    A = log((r1+r4)/(r2+r3))
    # entropy productions
    S = F*A
    return(S)
end

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    rates = [ k1*V, K1*X, k2*X*(X-1)*(X-2)/(V*V), K2*B*X*(X-1)/V]
    return(rates)
end

# function to calculate the time step
function timstep(rates::AbstractVector)
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# # function to advance gillespie one step
# function step(rates::AbstractVector,vars::Int64,reac::Int64)
#     r = rand()
#     rs = rates/sum(rates)
#     p = 0 # probability used for forward path
#     if r < rs[1]
#         vars += 1 # X produced from A (not simulated)
#         p = rs[1]
#         reac = 1
#     elseif r < rs[1] + rs[2]
#         vars -= 1 # X unravels to A (not simulated)
#         p = rs[2]
#         reac = 2
#     elseif r < rs[1] + rs[2] + rs[3]
#         vars -= 1 # X decays to B (not simulated)
#         p = rs[3]
#         reac = 3
#     else
#         vars += 1 # X regenerated from B (not simulated)
#         p = rs[4]
#         reac = 4
#     end
#     return(vars,p,reac)
# end
#
# # function to find reverse probability
# function rev(rs::AbstractVector,reac::Int64)
#     rs = rs/sum(rs)
#     if reac > 0 || reac < 5
#         if reac % 2 == 0 # even rates backwards
#             return(rs[reac-1])
#         else # odd forward
#             return(rs[reac+1])
#         end
#     else
#         error("Invalid reaction code returned")
#     end
# end

# function to advance gillespie one step
function step(rates::AbstractVector,vars::Int64,reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars += 1 # X produced from A (not simulated)
        p = rs[1] + rs[4]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels to A (not simulated)
        p = rs[2] + rs[3]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars -= 1 # X decays to B (not simulated)
        p = rs[3] + rs[2]
        reac = 3
    else
        vars += 1 # X regenerated from B (not simulated)
        p = rs[4] + rs[1]
        reac = 4
    end
    return(vars,p,reac)
end

# function to find reverse probability
function rev(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 5
        if reac == 2 || reac == 3 # even rates backwards
            return(rs[2]+rs[3])
        else # odd forward
            return(rs[1]+rs[4])
        end
    else
        error("Invalid reaction code returned")
    end
end

function Gillespie!(stead::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},times::Array{Float64,1},vars::Array{Int64,1},
                    V::Int64,reacs::Array{Int64,1})
    # Preallocate for simulation
    times[1] = 0
    vars[1] = stead
    reac = 0
    tup = tdown = 0
    for i = 1:noits
        # calculate rates
        rs = rates(vars[i],k1,K1,k2,K2,B,V)
        if i != 1
            # now use reac to calculate reverse rate
            pb[i-1] = rev(rs,reac)
        end
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[i+1], pf[i], reac = step(rs,vars[i],reac)
        posX = vars[i]# (vars[i+1] + vars[i])/(2)
        # final reverse rate
        if i == noits
            rs = rates(vars[end],k1,K1,k2,K2,B,V)
            pb[end] = rev(rs,reac)
        end
        reacs[i] = reac
    end
    return(pf,pb,times,vars,reacs)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,
                    star1::Float64,mid1::Float64,fin1::Float64,V::Int64)
    SX = zeros(noits)
    reacs = zeros(Int64,noits)
    # generate high and low states
    star2 = round(Int64,star1)
    fin2 = round(Int64,fin1)
    mid2 = round(Int64,mid1)
    # preallocating arrays used inside function
    pfX = zeros(noits)
    pbX = zeros(noits)
    timesX = zeros(noits+1)
    varsX = fill(0,noits+1)
    init = star2
    pfX, pbX, timesX, varsX, reacs = Gillespie!(init,k1,K1,k2,K2,B,noits,pfX,pbX,timesX,varsX,V,reacs)
    # calculate total entropy production
    println(minimum(varsX))
    println(maximum(varsX))
    for i = 1:noits
        SX[i] += log(pfX[i]) - log(pbX[i])
    end
    println("Gillespies Done!")
    return(SX,varsX,reacs,timesX[end])
end

# main function
function main()
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    B = 4.0
    V = 1500
    high2low = true
    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(k1,K1,k2,K2,B,high2low,V)
    # now run multiple Gillespie simulations
    noits = 25000000
    noruns = 1#500
    SX, varsX, reacs, t = multgill(noits,noruns,k1,K1,k2,K2,B,star1,mid1,fin1,V)
    # should do something cleverer
    Xmin = minimum(varsX)
    Xmax = maximum(varsX)
    L = Xmax - Xmin + 1
    S1 = zeros(L)
    S2 = zeros(L)
    S3 = zeros(L)
    S4 = zeros(L)
    Xrange = Xmin:Xmax
    for i = 1:length(SX)
        # Possibly should be varsX[i+1]
        if reacs[i] == 1
            S1[varsX[i]-Xmin+1] += SX[i]
        elseif reacs[i] == 2
            S2[varsX[i]-Xmin+1] += SX[i]
        elseif reacs[i] == 3
            S3[varsX[i]-Xmin+1] += SX[i]
        else
            S4[varsX[i]-Xmin+1] += SX[i]
        end
    end
    plot(Xrange,S1)
    savefig("../Results/test1.png")
    plot(Xrange,S2)
    savefig("../Results/test2.png")
    plot(Xrange,S3)
    savefig("../Results/test3.png")
    plot(Xrange,S4)
    savefig("../Results/test4.png")
    plot(Xrange,S1.+S4)
    savefig("../Results/test5.png")
    plot(Xrange,S2.+S3)
    savefig("../Results/test6.png")
    plot(Xrange,S1.+S2.+S3.+S4)
    savefig("../Results/test7.png")
    println(sum(SX)/(V*noruns*t))
    S = Sch(k1,K1,k2,K2,B,V,star1)
    Sr = Schr(k1,K1,k2,K2,B,V,star1)
    println(S/V)
    println(Sr/V)
    return(nothing)
end

@time main()
