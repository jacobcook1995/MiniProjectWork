#!/usr/bin/env julia
# gill2.jl
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
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1},Ω::Int64)
    # change this so that it uses less memory
    times = zeros(2)
    vars = fill(0,2,2)
    vars[:,2] = star
    hist = zeros(4*Ω,4*Ω)
    for i = 1:noits
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
        # add to histogram
        hist[vars[1,1]+1,vars[2,1]+1] += times[2] - times[1]
    end
    hist = hist/times[2]
    println("Gillespie Done!")
    return(hist)
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
    noits = 2500000000
    # maybe randomise the starting point somewhat
    # also maybe remove the vars
    hist = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,fin2,Ω)
    # now need to make matrices for W, rates of hopping between states
    # preallocating a sparse matrix to save space
    h1 = size(hist,1)
    h2 = size(hist,2)
    histw = spzeros(h1*h2,h1*h2)
    # now need to actually generate this matrix
    for j = 1:h2
        for i = 1:h1
            rs, dA, dB = rates(i-1,j-1,k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,true)
            histw[i+(j-1)*h1,i+(j-1)*h1] = -sum(rs) # self term flow out of state
            # flow to other states
            for l = 1:length(rs)
                # just cancel out flow from state if it's going out of considered region
                if dA[l] + i > h1 || dB[l] + j > h2
                    histw[i+(j-1)*h1,i+(j-1)*h1] += rs[l]
                # need to think about how to deal with negative cases
                elseif dA[l] + i < 1 || dB[l] + j < 1
                    if rs[l] != 0
                        error("Non-zero flux to negative part")
                    end
                else
                    histw[i+(j-1)*h2,i+dA[l]+(j+dB[l]-1)*h2] += rs[l]
                end
            end
        end
    end
    sm = spzeros(h1*h2,h1*h2)
    for j = 1:h2
        for i = 1:h1
            for k = -1:1
                if i + k <= h1 && i + k >= 1
                    sm[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+k+(j-1)*h2]*log(histw[i+(j-1)*h2,i+k+(j-1)*h2]/histw[i+k+(j-1)*h2,i+(j-1)*h2])
                end
                if j + k <= h2 && j + k >= 1
                    sm[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+(j+k-1)*h2]*log(histw[i+(j-1)*h2,i+(j+k-1)*h2]/histw[i+(j+k-1)*h2,i+(j-1)*h2])
                end
            end
        end
    end
    println(sum(sm))
    st = spzeros(h1*h2,h1*h2)
    for j = 1:h2
        for i = 1:h1
            for k = -1:1
                if i + k <= h1 && i + k >= 1 && hist[i,j] > 0 && hist[i+k,j] > 0
                    st[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+k+(j-1)*h2]*log(hist[i,j]*histw[i+(j-1)*h2,i+k+(j-1)*h2]/(hist[i+k,j]*histw[i+k+(j-1)*h2,i+(j-1)*h2]))
                end
                if j + k <= h2 && j + k >= 1 && hist[i,j] > 0 && hist[i,j+k] > 0
                    st[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+(j+k-1)*h2]*log(hist[i,j]*histw[i+(j-1)*h2,i+(j+k-1)*h2]/(hist[i,j+k]*histw[i+(j+k-1)*h2,i+(j-1)*h2]))
                end
            end
        end
    end
    println(sum(st))
    s = spzeros(h1*h2,h1*h2)
    for j = 1:h2
        for i = 1:h1
            for k = -1:1
                if i + k <= h1 && i + k >= 1 && hist[i,j] > 0 && hist[i+k,j] > 0
                    s[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+k+(j-1)*h2]*log(hist[i,j]/(hist[i+k,j]))
                end
                if j + k <= h2 && j + k >= 1 && hist[i,j] > 0 && hist[i,j+k] > 0
                    s[i,j] += hist[i,j]*histw[i+(j-1)*h2,i+(j+k-1)*h2]*log(hist[i,j]/(hist[i,j+k]))
                end
            end
        end
    end
    println(sum(s))
    sL = s[1:(mid2[1]+1),(mid2[2]+2):end]
    println(sum(sL))
    sH = s[(mid2[1]+2):end,1:(mid2[2]+1)]
    println(sum(sH))
    smL = sm[1:(mid2[1]+1),(mid2[2]+2):end]
    println(sum(smL))
    smH = sm[(mid2[1]+2):end,1:(mid2[2]+1)]
    println(sum(smH))
    stL = st[1:(mid2[1]+1),(mid2[2]+2):end]
    println(sum(stL))
    stH = st[(mid2[1]+2):end,1:(mid2[2]+1)]
    println(sum(stH))
    # now need to print histograms
    histA = zeros(h1) # each row reprents an A?
    for i = 1:h1
        histA[i] = sum(hist[i,:])
    end
    histB = zeros(h2) # each row reprents a B?
    for i = 1:h2
        histB[i] = sum(hist[:,i])
    end
    bar(histA)
    annotate!(2*Ω,0.666*maximum(histA),text("dS_t = $(sum(st))\ndS = $(sum(s))",:left))
    savefig("../Results/HistA$(ARGS[1]).png")
    bar(histB)
    annotate!(2*Ω,0.666*maximum(histB),text("dS_m = $(sum(sm))",:left))
    savefig("../Results/HistB$(ARGS[1]).png")
    # # new histograms
    # histL = hist[1:(mid2[1]+1),(mid2[2]+2):end]
    # histH = hist[(mid2[1]+2):end,1:(mid2[2]+1)]
    # # renormalise
    # histL = histL/sum(histL)
    # histH = histH/sum(histH)
    # SL = 0
    # for j = 1:size(histL,2)
    #     for i = 1:size(histL,1)
    #         if histL[i,j] != 0
    #             SL -= histL[i,j]*log(histL[i,j])
    #         end
    #     end
    # end
    # SH = 0
    # for j = 1:size(histH,2)
    #     for i = 1:size(histH,1)
    #         if histH[i,j] != 0
    #             SH -= histH[i,j]*log(histH[i,j])
    #         end
    #     end
    # end
    # S = 0
    # for j = 1:size(hist,2)
    #     for i = 1:size(hist,1)
    #         if hist[i,j] != 0
    #             S -= hist[i,j]*log(hist[i,j])
    #         end
    #     end
    # end
    # println(S)
    # println(Ω*SL)
    # println(Ω*SH)
    # println(Ω*(SH-SL))
    return(nothing)
end

@time main()
