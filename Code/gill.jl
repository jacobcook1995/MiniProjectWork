#!/usr/bin/env julia
# gill.jl
# A script to efficently perform a gillespie simulation of the full 4 species model
# This script generates this as a trajectory so that steady state entropy production
# can be calculated

using Roots
using Plots
import GR
# function to write out data as a csv
function writeout(times::Array{Float64,1},vars::Array{Int64,2})
    filename = "../Results/$(ARGS[1]).csv"
    out_file = open(filename, "w")
    # open file for writing
    for i = 1:size(vars,2)
        line = "$(vars[1,i]),$(vars[2,i]),$(vars[3,i]),$(vars[4,i]),$(times[i])\n"
        write(out_file, line)
    end
    # then close file
    close(out_file)
    return(nothing) # return nothing
end
# function to find the zeros of the function
function nullcline(F::Float64,r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,
                    q::Float64,kmin::Float64,qmin::Float64,high2low::Bool,Ne::Float64)
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
    sad = [ As[2]; Bs[2]; Ss[2]; Ws[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ss[1]; Ws[1] ]
        ss2 = [ As[3]; Bs[3]; Ss[3]; Ws[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ss[3]; Ws[3] ]
        ss2 = [ As[1]; Bs[1]; Ss[1]; Ws[1] ]
    end
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
function timstep(rates::AbstractVector)
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::AbstractVector,vars::Array{Int64,1})
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
                    f::Float64,r::Float64,F::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1})
    times = zeros(noits+1)
    vars = fill(0,4,noits+1)
    vars[:,1] = star
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1,i],vars[2,i],vars[3,i],vars[4,i],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,F)
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[:,i+1] = step(rs,vars[:,i])
    end
    println("Gillespie Done!")
    return(vars,times)
end

# main function
function main()
    # General parameters
    Ω = 60 # Ω = 1, gMAP parameterisation
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
    high2low = true

    # first need to use these parameters to find a steady state
    star1, _, _ = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low,Ne)
    # round star so that it becomes a vector of integers
    star2 = fill(0,4)
    for i = 1:4
        star2[i] = round(Int64,star1[i])
    end
    # now run gillespie
    noits = 5000000
    vars, times = gillespie(K,k,Q,q,kmin,qmin,f,r,F,Kmin,Qmin,noits,star2)
    if length(ARGS) >= 1
        writeout(times,vars)
    end

end

@time main()
