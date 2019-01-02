#!/usr/bin/env julia
# gill2spec.jl
# A script that performs a gillespie simulation of the 4 species genetic switch model
# with explicit promotor un/binding. Histograms of the occupation and the waiting time
# distributions of the two states are generated and output
#
# Author: Jacob Cook
# Date: September 2018

using Roots
using SymPy

# function to provide normalised parameters for the gillespie simulation
function paras(Ω::Int64)
    # just gonna hard code some in for the moment
    # EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11.0/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/(Ω^2) # Promoter switching
    r = 10.0
    F = 10.0*Ω
    Kmin = 10.0^-10 # remains neligable though
    Qmin = 10.0^-10
    Ne = 150.0*Ω # number of elements in the system
    return(k,kmin,q,qmin,K,Kmin,Q,Qmin,r,f,F,Ne)
end

# ps = [K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,Ne]
# function to generate steady states
function states(ps::Array{Float64,1},Ω::Int64,maxr::Float64,δ::Float64)
    # unnormalise the parameters as this is at Langevin level
    ps[12] = ps[12]/Ω
    ps[11] = ps[11]/Ω
    ps[9] = ps[9]*(Ω^2)
    # define symbols
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N = SymPy.symbols("A,B,S,W,K k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,N")
    # equation for S
    Se = N - W - A - B
    # equation for W
    We = (K*A + Q*B - F)/(Qmin + Kmin)
    # eliminate W from equation for S
    Se = SymPy.subs(Se,W=>We)
    # equation for dB/dt
    dB = (q*S*r)/(r + f*A^2) - qmin*B - Q*B + Qmin*W
    # sub in S and W expressions
    dB = SymPy.subs(dB,S=>Se,W=>We)
    # solve this expression for B
    Bear = solve(dB,B)
    # remove from being in array form
    Be = Bear[1]
    # now sub this into Se and We
    Se = SymPy.subs(Se,B=>Be)
    We = SymPy.subs(We,B=>Be)
    # equation for dA/dt
    dA = (k*S*r)/(r + f*B^2) - kmin*A - K*A + Kmin*W
    # sub in B, S and W expressions
    dA = SymPy.subs(dA,B=>Be,S=>Se,W=>We)
    # then sub in parameters
    dA = SymPy.subs(dA,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    dA = SymPy.subs(dA,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # now find three stationary points in A
    guess = 0
    three = false
    As = Array{Float64,1}(undef,0)
    n = 0
    while three == false
        a = convert(Float64,nsolve(dA,guess))
        # add to vector if not already found
        if (a in As) == false
            As = vcat(As,a)
            n += 1
        end
        if n == 3
            three = true
        end
        guess += δ
        if guess >= maxr
            println("Could not find three stationary points in range (0,$(guess)) with step $(δ)")
            println(As)
            error()
        end
    end
    # sub parameters into Be, Se, We
    Be = SymPy.subs(Be,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Be = SymPy.subs(Be,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    Se = SymPy.subs(Se,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Se = SymPy.subs(Se,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    We = SymPy.subs(We,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    We = SymPy.subs(We,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # then sub the various As values to get final values
    Bs = zeros(3)
    Ss = zeros(3)
    Ws = zeros(3)
    for i = 1:3
        Bs[i] = SymPy.subs(Be,A=>As[i]) |> float
        Ss[i] = SymPy.subs(Se,A=>As[i]) |> float
        Ws[i] = SymPy.subs(We,A=>As[i]) |> float
    end
    # finally put into states
    ss1 = [As[1],Bs[1],Ss[1],Ws[1]]
    sad = [As[2],Bs[2],Ss[2],Ws[2]]
    ss2 = [As[3],Bs[3],Ss[3],Ws[3]]
    println(ss1)
    println(sad)
    println(ss2)
    flush(stdout)
    return(ss1,sad,ss2)
end

# function to construct the rates for non-explicit switching case
function rates(A::Int64,B::Int64,S::Int64,W::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64,
                F::Float64,diffs::Bool=false)
    rates = [ k*S*r/(r+f*B^2), kmin*A, K*A, Kmin*W, q*S*r/(r+f*A^2), qmin*B, Q*B, Qmin*W, F]
    if diffs == false
        return(rates)
    else
        dA = [ 1, -1, -1, 1, 0, 0, 0, 0, 0 ]
        dB = [ 0, 0, 0, 0, 1, -1, -1, 1, 0 ]
        dS = [ -1, 1, 0, 0, -1, 1, 0, 0, 1 ]
        dW = [ 0, 0, 1, -1, 0, 0, 1, -1, -1 ]
        return(rates,dA,dB,dS,dW)
    end
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step: explicit switching
function step!(rates::Array{Float64,1},vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars[1] += 1 # A produced
        vars[3] -= 1 # S used
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
        vars[3] += 1 # S regenerated
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
        vars[4] += 1 # W formed
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1 # A regenerated
        vars[4] -= 1 # W removed
    elseif r < sum(rs[1:5])
        vars[2] += 1 # B produced
        vars[3] -= 1 # S used
    elseif r < sum(rs[1:6])
        vars[2] -= 1 # B unravels
        vars[3] += 1 # S regenerated
    elseif r < sum(rs[1:7])
        vars[2] -= 1 # B decays
        vars[4] += 1 # W formed
    elseif r < sum(rs[1:8])
        vars[2] += 1 # B regenerated
        vars[4] -= 1 # W removed
    else
        vars[3] += 1 # S supplied
        vars[4] -= 1 # W removed
    end
    return(vars) # there's always the potential of inverting this to speed the code up
end

# function to actually run a gillepsie simulation
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,F::Float64,Ne::Float64,
                    noits::Int64,star::Array{Int64,1},Ω::Int64,fin::Array{Int64,1})
    # define high A and high B states for comparison
    hA = fin[1]
    hB = star[2]
    # set bool as being in high B state initially
    highA = false
    # make vector to store waiting times
    wA = Array{Float64,1}(undef,0)
    wB = Array{Float64,1}(undef,0)
    tp = 0 # prior time
    times = zeros(2)
    vars = fill(0,4,2)
    vars[:,2] = star
    vol = convert(Int64,Ne+1)
    hist = zeros(35*Ω,35*Ω)
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
        vars[:,2] = step!(rs,vars[:,1])
        # add to histogram
        hist[vars[1,1]+1,vars[2,1]+1] += times[2] - times[1]
        # now step to check if transistion has occurred
        if vars[1,2] >= hA && highA == false
            highA = true
            wB = vcat(wB,times[2]-tp)
            tp = times[2]
        elseif vars[2,2] >= hB && highA == true # inefficent
            highA = false
            wA = vcat(wA,times[2]-tp)
            tp = times[2]
        end
    end
    hist = hist/times[2]
    return(hist,wA,wB,times[2])
end

# function to printout my data
function printout(hist::Array{Float64,2},wA::Array{Float64,1},wB::Array{Float64,1},Time::Float64)
    if isfile("../Results/wA$(ARGS[1])V$(ARGS[2]).csv")
        # if file already exists should add to bottom of it
        output_file = "../Results/wA$(ARGS[1])V$(ARGS[2]).csv"
        out_file = open(output_file, "a")
        for i = 1:size(wA,1)
            line = "$(wA[i])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        output_file = "../Results/wA$(ARGS[1])V$(ARGS[2]).csv"
        out_file = open(output_file, "w")
        for i = 1:size(wA,1)
            line = "$(wA[i])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    if isfile("../Results/wB$(ARGS[1])V$(ARGS[2]).csv")
        # if file already exists should add to bottom of it
        output_file = "../Results/wB$(ARGS[1])V$(ARGS[2]).csv"
        out_file = open(output_file, "a")
        for i = 1:size(wB,1)
            line = "$(wB[i])\n"
            write(out_file,line)
        end
        close(out_file)
    else
        # if file doesn't exist should make it
        output_file = "../Results/wB$(ARGS[1])V$(ARGS[2]).csv"
        out_file = open(output_file, "w")
        for i = 1:size(wB,1)
            line = "$(wB[i])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    # histogram in this case is two dimensional so not hard to represent
    if isfile("../Results/hist$(ARGS[1])V$(ARGS[2]).csv")
        # need to read in old file and then add new histogram to it and then output
        input_file = "../Results/hist$(ARGS[1])V$(ARGS[2]).csv"
        vol = convert(Int64,(countlines(input_file)) - 3)
        histr = Array{Float64,2}(undef,vol,vol)
        a = 0
        i = 1
        j = 1
        T = 0
        open(input_file, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                if line[1] == '#'
                    a += 1
                    i = 1
                    j = 1
                else
                    if a == 2
                        T = parse(Float64,line)
                    else
                        # need to think about how to do the parsing here
                        L = length(line)
                        comma = fill(0,vol+1)
                        j = 1
                        for i = 1:L
                            if line[i] == ','
                                j += 1
                                comma[j] = i
                            end
                        end
                        comma[end] = L+1
                        for j = 1:vol
                            histr[i,j] = parse(Float64,line[(comma[j]+1):(comma[j+1]-1)])
                        end
                        i += 1
                    end
                end
            end
        end
        hist = (Time*hist .+ T*histr)/(T + Time)
        # update time to include time of new run
        Time += T
    end
    # if file doesn't exist should make it
    output_file = "../Results/hist$(ARGS[1])V$(ARGS[2]).csv"
    out_file = open(output_file, "w")
    line = "# Start #\n"
    write(out_file,line)
    for i = 1:size(hist,1)
        line = ""
        for j = 1:size(hist,2)
            line *= "$(hist[i,j,1,1]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    line = "# Time #\n"
    write(out_file,line)
    line = "$(Time)\n"
    write(out_file,line)
    close(out_file)
    return(nothing)
end

function main()
    println("Compiled, Starting script.")
    flush(stdout)
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
    k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, F, Ne = paras(Ω)
    # use these parameters to generate the steady states
    ss1, sad, ss2 = states([K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,Ne],Ω,15.0,0.1)
    # scale starting posisition by volume
    star = fill(0,4)
    fin = fill(0,4)
    for i = 1:4
        star[i] = round(Int64,ss1[i]*Ω)
        fin[i] = round(Int64,ss2[i]*Ω)
    end
    # Now ready to set the gillespie simulation running
    noits = 10000000000
    hist, wA, wB, T = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,F,Ne,noits,star,Ω,fin)
    # now printout
    printout(hist,wA,wB,T)
    ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, F, Ne, Ω ]
    output_file = "../Results/ps$(ARGS[1])V$(ARGS[2]).csv"
    out_file = open(output_file, "w")
    for i = 1:length(ps)
        line = "$(ps[i])\n"
        write(out_file,line)
    end
    close(out_file)
    return(nothing)
end

@time main()
