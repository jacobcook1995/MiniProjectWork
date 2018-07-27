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
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,F::Float64,diffs::Bool=false)
    rates = [r*k*S/(r+f*B*(B-1)), kmin*A, K*A, Kmin*W, r*q*S/(r+f*A*(A-1)), qmin*B, Q*B, Qmin*W, F]
    if diffs == false
        return(rates)
    else
        dA = [ 1, -1, -1, 1, 0, 0, 0, 0, 0]
        dB = [ 0, 0, 0, 0, 1, -1, -1, 1, 0]
        dS = [ -1, 1, 0, 0, -1, 1, 0, 0, 1]
        dW = [ 0, 0, 1, -1, 0, 0, 1, -1, -1]
        return(rates,dA,dB,dS,dW)
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
    hist = zeros(3*maxA,3*maxB,3*maxS,3*maxW)
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
        hist[vars[1,1]+1,vars[2,1]+1,vars[3,1]+1,vars[4,1]+1] += times[2] - times[1]
    end
    hist = hist/times[2]
    println("Gillespie Done!")
    return(hist)
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
    noits = 50#000000
    # maybe randomise the starting point somewhat
    # also maybe remove the vars
    hist = gillespie(K,k,Q,q,kmin,qmin,f,F,r,Kmin,Qmin,noits,star2,maxA,maxB,maxS,maxW)
    h1 = size(hist,1)
    h2 = size(hist,2)
    h3 = size(hist,3)
    h4 = size(hist,4)
    # preallocate as a sparse array
    h12 = h1*h2
    h123 = h12*h3
    h1234 = h123*h4
    @time I = zeros(Int64,10*h1234)
    j = 1
    for i = 1:(10*h1234)
        I[i] .= j
        if i % 10 == 0
            j += 1
        end
    end
    @time J = zeros(10*h1234)
    tic()
    for i = 1:h1234
        # first decompose i, this is wrong
        w = floor((i-1)/h123)
        s = floor((i-1-w*h123)/h12)
        b = floor((i-1-w*h123-s*h12)/h1)
        a = i-b*h1-s*h12-w*h123-1
        J[1+(i-1)*10] .= i
        # check a and s in valid range for reac 1
        if a == 0 || s >= h3
            J[3+(i-1)*10] .= 0
        else
            J[3+(i-1)*10] .= i - 1 + h12
        end
        if a == 0 || w >= h4
            J[6+(i-1)*10] .= 0
        else
            J[6+(i-1)*10] .= i - 1 + h123
        end
        if a >= h1 || s == 0
            J[2+(i-1)*10] .= 0
        else
            J[2+(i-1)*10] .= i + 1 - h12
        end
        if a >= h1 || w == 0
            J[7+(i-1)*10] .= 0
        else
            J[7+(i-1)*10] .= i + 1 - h123
        end
        if b == 0 || s >= h3
            J[5+(i-1)*10] .= 0
        else
            J[5+(i-1)*10] .= i - h1 + h12
        end
        if b == 0 || w >= h4
            J[8+(i-1)*10] .= 0
        else
            J[8+(i-1)*10] .= i - h1 + h123
        end
        if b >= h2 || s == 0
            J[4+(i-1)*10] .= 0
        else
            J[4+(i-1)*10] .= i + h1 - h12
        end
        if b >= h2 || w == 0
            J[9+(i-1)*10] .= 0
        else
            J[9+(i-1)*10] .= i + h1 - h123
        end
        # deal with supply case
        if s >= h3 || w == 0
            J[10+(i-1)*10] .= 0
        else
            J[10+(i-1)*10] .= i + h12 - h123
        end
    end
    toc()
    @time t = find(J .== 0)
    # now need to remove zeros elements from J and corresponding elements from I
    J = deleteat!(J, t)
    I = deleteat!(I, t)
    @time V = zeros(length(J))
    gc()
    histw = sparse(I,J,V)
    return(nothing)
    # now need to actually generate this matrix
    # This is clearly not going to work as the matrix will be too large to actually store
    for m = 1:h4
        for l = 1:h3
            for j = 1:h2
                for i = 1:h1
                    rs, dA, dB, dS, dW = rates(i-1,j-1,l-1,m-1,k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f,F,true)
                    pos = i + (j-1)*h1 + (l-1)*h12 + (m-1)*h123
                    histw[pos,pos] .= sum(rs) # self term flow out of state
                    # flow to other states
                    for n = 1:length(rs)
                        # just cancel out flow from state if it's going out of considered region
                        if dA[n] + i > h1 || dB[n] + j > h2 || dS[n] + l > h3 || dW[n] + m > h4
                            histw[pos,pos] += rs[n] #@time
                        elseif dA[n] + i < 1 || dB[n] + j < 1 || dS[n] + l < 1 || dW[n] + m < 1
                            if m == 1 && n == 9
                                histw[pos,pos] += rs[n]
                            elseif l == h3 && n == 9
                                histw[pos,pos] += rs[n]
                            elseif rs[n] != 0
                                println(i,j,l,m,n)
                                println(dA[n],dB[n],dS[n],dW[n])
                                error("Non-zero flux to negative part")
                            end
                        else
                            histw[pos,pos + dA[n]+ dB[n]*h1 + dS[n]*h12 + dW[n]*h123] += rs[n] #@time
                        end

                    end
                end
            end
        end
    end
    # now split into high and low histograms
    histL = hist[1:(mid2[1]+1),(mid2[2]+2):end,:,:]
    histH = hist[(mid2[1]+2):end,1:(mid2[2]+1),:,:]
    hist1 = hist[1:(mid2[1]+1),1:(mid2[2]+1),:,:]
    hist2 = hist[(mid2[1]+2):end,(mid2[2]+2):end,:,:]
    # renormalise
    # histL = histL/sum(histL)
    # histH = histH/sum(histH)
    # calculate entropy
    SL = 0
    for m = 1:size(histL,4)
        for l = 1:size(histL,3)
            for j = 1:size(histL,2)
                for i = 1:size(histL,1)
                    if histL[i,j,l,m] != 0
                        SL -= histL[i,j,l,m]*log(histL[i,j,l,m])
                    end
                end
            end
        end
    end
    SH = 0
    for m = 1:size(histH,4)
        for l = 1:size(histH,3)
            for j = 1:size(histH,2)
                for i = 1:size(histH,1)
                    if histH[i,j,l,m] != 0
                        SH -= histH[i,j,l,m]*log(histH[i,j,l,m])
                    end
                end
            end
        end
    end
    S = 0
    for m = 1:size(hist,4)
        for l = 1:size(hist,3)
            for j = 1:size(hist,2)
                for i = 1:size(hist,1)
                    if hist[i,j,l,m] != 0
                        S -= hist[i,j,l,m]*log(hist[i,j,l,m])
                    end
                end
            end
        end
    end
    S1 = 0
    for m = 1:size(hist1,4)
        for l = 1:size(hist1,3)
            for j = 1:size(hist1,2)
                for i = 1:size(hist1,1)
                    if hist1[i,j,l,m] != 0
                        S1 -= hist1[i,j,l,m]*log(hist1[i,j,l,m])
                    end
                end
            end
        end
    end
    S2 = 0
    for m = 1:size(hist2,4)
        for l = 1:size(hist2,3)
            for j = 1:size(hist2,2)
                for i = 1:size(hist2,1)
                    if hist2[i,j,l,m] != 0
                        S2 -= hist2[i,j,l,m]*log(hist2[i,j,l,m])
                    end
                end
            end
        end
    end
    println(S)
    println(SL)
    println(SH)
    println(S-SL-SH)
    println(S1)
    println(S2)
    println(S1+S2)
    return(nothing)
end

@time main()
