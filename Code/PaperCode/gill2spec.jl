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
    flush(stdout) # needed to ensure output prints when using nohup
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
        vars[2] += 2 # relasing 2 B
    elseif r < sum(rs[1:10])
        vars[4] += 1 # promotor b unbinds
        vars[1] += 2 # releasing two A
    elseif r < sum(rs[1:11])
        vars[3] -= 1 # promotor a binds
        vars[2] -= 2 # bound by 2 B
    else
        vars[4] -= 1 # promotor b binds
        vars[1] -= 2 # bound by 2 A
    end
    return(vars) # there's always the potential of inverting this to speed the code up
end

# function to actually run a gillepsie simulation
function gillespie(K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,qmin::Float64,
                    f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,star::Array{Int64,1},
                    Ω::Int64,fin::Array{Int64,1})
    # define high A and high B states for comparison
    hA = star[1]
    hB = fin[2]
    # set bool as being in high A state initially
    highA = true
    # make vector to store waiting times
    wA = Array{Float64,1}(undef,0)
    wB = Array{Float64,1}(undef,0)
    tp = 0 # prior time
    times = zeros(2)
    vars = fill(0,4,2)
    vars[:,2] = star
    hist = zeros(25*Ω,25*Ω,2,2)
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
function printout(hist::Array{Float64,4},wA::Array{Float64,1},wB::Array{Float64,1},Time::Float64)
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
    # histogram is more tricky as have to represent 4D array as a 2D text file
    if isfile("../Results/hist$(ARGS[1])V$(ARGS[2]).csv")
        # need to read in old file and then add new histogram to it and then output
        input_file = "../Results/hist$(ARGS[1])V$(ARGS[2]).csv"
        vol = convert(Int64,(countlines(input_file)/4) - 1.5)
        histr = Array{Float64,4}(undef,vol,vol,2,2)
        a = 0
        i = 1
        j = 1
        k = 1
        l = 1
        T = 0
        open(input_file, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                if line[1] == '#'
                    a += 1
                    i = 1
                    j = 1
                    if l == 1 && a != 1
                        l = 2
                    elseif a != 1
                        l = 1
                        k = 2
                    end
                else
                    if a == 5
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
                            histr[i,j,k,l] = parse(Float64,line[(comma[j]+1):(comma[j+1]-1)])
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
    line = "# (0-0) #\n"
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
    line = "# (0-1) #\n"
    write(out_file,line)
    for i = 1:size(hist,1)
        line = ""
        for j = 1:size(hist,2)
            line *= "$(hist[i,j,1,2]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    line = "# (1-0) #\n"
    write(out_file,line)
    for i = 1:size(hist,1)
        line = ""
        for j = 1:size(hist,2)
            line *= "$(hist[i,j,2,1]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    line = "# (1-1) #\n"
    write(out_file,line)
    for i = 1:size(hist,1)
        line = ""
        for j = 1:size(hist,2)
            line *= "$(hist[i,j,2,2]),"
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
    noits = 10000000000
    hist, wA, wB, T = gillespie(K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,star,Ω,fin)
    # now printout
    printout(hist,wA,wB,T)
    ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f, Ω ]
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
