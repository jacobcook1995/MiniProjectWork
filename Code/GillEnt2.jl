#!/usr/bin/env julia
# GillEnt2.jl
# A script to perform a Gillespie simulation and then calculate the entropy production
# of the paths by direct computation of the reversed probabilities.
# The script will then calculate steady state entropy via the Shannon method for comparison
# This will first be done for the two species case

using Roots

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
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to calculate the time step
function timstep(rates::AbstractVector)
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::AbstractVector,vars::Array{Int64,1},reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars[1] += 1 # A produced
        p = rs[1]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
        p = rs[2]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
        p = rs[3]
        reac = 3
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1 # A regenerated
        p = rs[4]
        reac = 4
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5]
        vars[2] += 1 # B produced
        p = rs[5]
        reac = 5
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6]
        vars[2] -= 1 # B unravels
        p = rs[6]
        reac = 6
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7]
        vars[2] -= 1 # B decays
        p = rs[7]
        reac = 7
    else
        vars[2] += 1 # B regenerated
        p = rs[8]
        reac = 8
    end
    return(vars,p,reac)
end

# function to find reverse probability
function rev(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 9
        if reac % 2 == 0 # even rates backwards
            return(rs[reac-1])
        else # odd forward
            return(rs[reac+1])
        end
    else
        error("Invalid reaction code returned")
    end
end

function Gillespie!(stead::Array{Int64,1},K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,
                    qmin::Float64,f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},times::Array{Float64,1},vars::Array{Int64,2})
    # Preallocate for simulation
    times[1] = 0
    vars[:,1] = stead
    reac = 0
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1,i],vars[2,i],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        if i != 1
            # now use reac to calculate reverse rate
            pb[i-1] = rev(rs,reac)
        end
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[:,i+1], pf[i], reac = step(rs,vars[:,i],reac)
        # final reverse rate
        if i == noits
            rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            pb[end] = rev(rs,reac)
        end
    end
    return(pf,pb,times,vars)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,r::Float64,f::Float64,K::Float64,Q::Float64,
                    k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64,
                    highA::Array{Int64,1},highB::Array{Int64,1})
    SA = zeros(noruns)
    SB = zeros(noruns)
    minA = zeros(noruns)
    minB = zeros(noruns)
    maxA = zeros(noruns)
    maxB = zeros(noruns)
    # preallocating arrays used inside function
    pfA = zeros(noits)
    pbA = zeros(noits)
    timesA = zeros(noits+1)
    varsA = fill(0,2,noits+1)
    pfB = zeros(noits)
    pbB = zeros(noits)
    timesB = zeros(noits+1)
    varsB = fill(0,2,noits+1)
    for j = 1:noruns
        pfA, pbA, timesA, varsA = Gillespie!(highA,K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,pfA,pbA,timesA,varsA)
        pfB, pbB, timesB, varsB = Gillespie!(highB,K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,pfB,pbB,timesB,varsB)
        # calculate total entropy production
        for i = 1:noits
            SA[j] += log(pfA[i]) - log(pbA[i])
            SB[j] += log(pfB[i]) - log(pbB[i])
        end

        # convert total entropy production to entropy production rate
        SA[j] = SA[j]/timesA[end]
        SB[j] = SB[j]/timesB[end]
        # store this data so validity can be checked later on
        minA[j] = minimum(varsA[1,:])
        minB[j] = minimum(varsB[2,:])
        maxA[j] = maximum(varsA[2,:])
        maxB[j] = maximum(varsB[1,:])
    end
    println("Gillespies Done!")
    return(SA,SB,minA,minB,maxA,maxB)
end

# function to compute the shannon entropy at a fixed point
function shannon(point::Array{Float64,1},r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,
                kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    pa = k*r/(r+f*(point[2]-1)*point[2])
    pamin = kmin*point[1]
    pb = q*r/(r+f*(point[1]-1)*point[1])
    pbmin = qmin*point[2]
    da = K*point[1]
    damin = Kmin
    db = Q*point[2]
    dbmin = Qmin
    S1 = (pa - pamin)*log(pa/pamin)
    S2 = (pb - pbmin)*log(pb/pbmin)
    S3 = (da - damin)*log(da/damin)
    S4 = (db - dbmin)*log(db/dbmin)
    S = S1 + S2 + S3 + S4
    return(S)
end

# main function
function main()
    # General parameters
    Ω = 300.0
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
    star1, _, fin1 = nullcline(r,f,K,Q,k,q,kmin,qmin)
    # round star so that it becomes a vector of integers
    star2 = fill(0,2)
    fin2 = fill(0,2)
    for i = 1:2
        star2[i] = round(Int64,star1[i])
        fin2[i] = round(Int64,fin1[i])
    end
    # now the steady state will be used to find shannon entropy productions
    SAS = shannon(star1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    SBS = shannon(fin1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println("Entropy production rate of high A state via Shannon formula = $(SAS)")
    println("Entropy production rate of high B state via Shannon formula = $(SBS)")
    # now run multiple Gillespie simulations
    noits = 500000
    noruns = 50
    SA, SB, minA, minB, maxA, maxB = multgill(noits,noruns,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin,star2,fin2)
    println(sum(SA)/noruns)
    println(sum(SB)/noruns)
    return(nothing)
end

@time main()
