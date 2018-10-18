#!/usr/bin/env julia
# gill2pl.jl
# A script to make nice plots of a two species gillespie simulation

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
    println(ss1)
    println(sad)
    println(ss2)
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
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1},Ω::Int64)
    # change this so that it uses less memory
    ts = zeros(noits+1)
    ABs = zeros(noits+1,2)
    t = 0
    AB = copy(star)
    ts[1] = t
    ABs[1,:] = AB
    for i = 1:noits
        # calculate rates
        rs = rates(AB[1],AB[2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # calculate timestep
        τ = timstep(rs)
        # update time
        t += τ
        # do gillepsie step
        AB = step(rs,AB)
        # add to vectors
        ts[i+1] = t
        ABs[i+1,:] = AB
    end
    println("Gillespie Done!")
    return(ABs,ts)
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
    r = 100.0

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
    noits = 250000000
    # maybe randomise the starting point somewhat
    # also maybe remove the vars
    ABs, ts = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,fin2,Ω)
    # should actually find the crossing points then just plot
    trig = false
    pos1 = 0
    while trig == false
        pos1 += 1
        if ABs[pos1,1] >= star2[1] && ABs[pos1,2] <= star2[2]
            trig = true
        elseif pos1 == noits
            error("Not long enough simulation\n")
        end
    end
    intv = 50000
    ABsr = ABs[(pos1-intv):(pos1),:]
    tsr = ts[(pos1-intv):(pos1)]
    trig = false
    pos2 = 0
    while trig == false
        pos2 += 1
        if ABsr[intv-pos2,2] >= fin2[2] && ABsr[intv-pos2,1] <= fin2[1]
            trig = true
        end
    end
    trig = false
    pos3 = pos1
    while trig == false
        pos3 += 1
        if ABs[pos3,2] >= fin2[2] && ABs[pos3,1] <= fin2[1]
            trig = true
        elseif pos3 == noits
            error("Not long enough simulation\n")
        end
    end
    ABsf = ABs[(pos3-intv):(pos3),:]
    tsf = ts[(pos3-intv):(pos3)]
    trig = false
    pos4 = 0
    while trig == false
        pos4 += 1
        if ABsf[intv-pos4,2] <= star2[2] && ABsf[intv-pos4,1] >= star2[1]
            trig = true
        end
    end
    plot(ABsr[(intv-pos2):end,2],ABsr[(intv-pos2):end,1],label="to high A",title="2D Toggle Switch")
    plot!(ABsf[(intv-pos4):end,2],ABsf[(intv-pos4):end,1],label="to high B",xlabel="A",ylabel="B")
    scatter!([star2[2]], [star2[1]], seriescolor = :green, label="")
    scatter!([fin2[2]], [fin2[1]], seriescolor = :red, label="")
    savefig("../Results/switch.png")
    return(nothing)
end

@time main()
