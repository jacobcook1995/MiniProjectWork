#!/usr/bin/env julia
# gill2spec.jl
# A script that performs a gillespie simulation of the 2 species genetic switch model
# with explicit promotor un/binding. Histograms of the occupation and the waiting time
# distributions of the two states are generated and output
#
# Author: Jacob Cook
# Date: September 2018

using Roots
using Plots
import GR # this is necessary to avoid a world age error when using GR in function

# function to provide normalised parameters for the gillespie simulation
function paras(Ω::Int64)
    # just gonna hard code some in for the moment
    # EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT
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
    return(k,kmin,q,qmin,K,Kmin,Q,Qmin,r,f)
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
    println("The three stationary points are:")
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end


# function to construct the rates
function rates(A::Int64,B::Int64,a::Int64,b::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,diffs::Bool=false)
    rates = [ k*a, kmin*A, K*A, Kmin, q*b, qmin*B, Q*B, Qmin, r*(1-a), r*(1-b), f*a*(B-1)*B, f*b*(A-1)*A ]
    if diffs == false
        return(rates)
    else
        dA = [ 1, -1, -1, 1, 0, 0, 0, 0, 0, 2, 0, -2 ]
        dB = [ 0, 0, 0, 0, 1, -1, -1, 1, 2, 0, -2, 0 ]
        da = [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0]
        db = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1]
        return(rates,dA,dB,da,db)
    end
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step!(rates::Array{Float64,1},vars::Array{Int64,1})
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
    elseif r < sum(rs[1:5])
        vars[2] += 1 # B produced
    elseif r < sum(rs[1:6])
        vars[2] -= 1 # B unravels
    elseif r < sum(rs[1:7])
        vars[2] -= 1 # B decays
    elseif r < sum(rs[1:8])
        vars[2] += 1 # B regenerated
    elseif r < sum(rs[1:9])
        vars[3] += 1 # promotor a unbinds
    elseif r < sum(rs[1:10])
        vars[4] += 1 # promotor b unbinds
    elseif r < sum(rs[1:11])
        vars[3] -= 1 # promotor a binds
    else
        vars[4] -= 1 # promotor b binds
    end
    return(vars) # there's always the potential of inverting this to speed the code up
end

# function to actually run a gillepsie simulation
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1},Ω::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,4,2)
    vars[:,2] = star
    hist = zeros(40*Ω,40*Ω,2,2)#(20*Ω,20*Ω)
    for i = 1:noits
        vars[:,1] = vars[:,2]
        times[1] = times[2]
        # calculate rates
        rs = rates(vars[1,1],vars[2,1],vars[3,1],vars[4,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[2] = times[1] + τ
        # do gillepsie step
        vars[:,2] = step!(rs,vars[:,1])
        # add to histogram
        hist[vars[1,1]+1,vars[2,1]+1,vars[3,1]+1,vars[4,1]+1] += times[2] - times[1]
    end
    hist = hist/times[2]
    return(hist)
end

function main()
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    # then check that a system volume has been provided
    elseif length(ARGS) == 1
        println("Error: Need to provide an argument to set system volume.")
        return(nothing)
    end
    # Then take system volume Ω, check if provided value is integer
    Ω = 0
    try Ω = parse(Int64,ARGS[2])
    catch y
        if isa(y, ArgumentError) # would only really expect an argument error
            println("Error: System volume has to be integer.")
            return(nothing)
        end
    end
    # now want to obtain the parameters of the simulation
    k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f = paras(Ω)
    # use these parameters to generate the steady states
    ss1, sad, ss2 = states(k,kmin,q,qmin,K,Kmin,Q,Qmin,r,f,Ω)
    # scale starting posisition by volume
    star = [ 0, 0, 1, 0 ]
    fin = [ 0, 0, 0, 1 ]
    for i = 1:2
        star[i] = round(Int64,ss1[i]*Ω)
        fin[i] = round(Int64,ss2[i]*Ω)
    end
    # Now ready to set the gillespie simulation running
    noits = 100000000
    hist = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,star,Ω)
    hist = dropdims(sum(hist,dims=4),dims=4)
    hist = dropdims(sum(hist,dims=3),dims=3)
    histA = dropdims(sum(hist,dims=2),dims=2)
    bar(histA)
    savefig("../../Results/SmartName.png")
    histB = dropdims(sum(hist,dims=1),dims=1)
    bar(histB)
    savefig("../../Results/SmartName2.png")
    return(nothing)
end

@time main()
