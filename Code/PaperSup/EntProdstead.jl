#!/usr/bin/env julia
# EntProdtraj.jl
# A script to read in parameters and steady states and then to generate nice histograms
# of the entropy productions by three differenet methods.
# 1) using the langevin entropy production term
# 2) by using a large volume reduced form of the master equation
# 3) by using the underlying master equation
#
# Author: Jacob Cook
# Date: January 2019

using Plots
using LaTeXStrings
using PyCall
import PyPlot

# Diffusion matrix D
# ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
function D!(D::Array{Float64,2},x::Array{Float64,1},ps::Array{Float64,1})
    D[1,1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) + (ps[5]+ps[2])*x[1] + ps[6]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) + (ps[7]+ps[4])*x[2] + ps[8]
    return(D)
end

# vector of forces
function f!(f::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    Ω = 5000
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) - ps[2]*x[1] - ps[5]*x[1] + ps[6]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) - ps[4]*x[2] - ps[7]*x[2] + ps[8]
    return(f)
end

# A function to find and plot the langevin entropy productions of the various trajectories
function LangEnt(traj::Array{Float64,2},ps::Array{Float64,1},ts::Array{Float64,1})
    entp = zeros(size(traj,1)-1)
    f = [0.0; 0.0]
    D = [0.0 0.0; 0.0 0.0]
    for i = 1:length(entp)
        dt = ts[i+1] - ts[i]
        qdot = (traj[i+1,:] .- traj[i,:])/(dt)
        # and need to find midpoint
        posA = (traj[i+1,1] + traj[i,1])/(2)
        posB = (traj[i+1,2] + traj[i,2])/(2)
        f = f!(f,[posA,posB],ps)
        D = D!(D,[posA,posB],ps)
        for j = 1:2
            for k = 1:2
                if D[j,k] != 0
                    entp[i] += 2*qdot[j]*f[k]*dt/D[j,k]
                end
            end
        end
    end
    entP = sum(entp)
    return(entP)
end

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to calculate the time step
function timstep(rates::Array{Float64,1})
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Array{Int64,1},reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars[1] += 1 # A produced
        p = rs[1]
        pt = rs[1] + rs[4]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
        p = rs[2]
        pt = rs[2] + rs[3]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
        p = rs[3]
        pt = rs[3] + rs[2]
        reac = 3
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1 # A regenerated
        p = rs[4]
        pt = rs[4] + rs[1]
        reac = 4
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5]
        vars[2] += 1 # B produced
        p = rs[5]
        pt = rs[5] + rs[8]
        reac = 5
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6]
        vars[2] -= 1 # B unravels
        p = rs[6]
        pt = rs[6] + rs[7]
        reac = 6
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7]
        vars[2] -= 1 # B decays
        p = rs[7]
        pt = rs[7] + rs[6]
        reac = 7
    else
        vars[2] += 1 # B regenerated
        p = rs[8]
        pt = rs[8] + rs[5]
        reac = 8
    end
    return(vars,p,pt,reac)
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

# function to find reverse probability
function rev(rs::Array{Float64,1},reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 9
        if reac == 1 # A produced
            return(rs[2],rs[2]+rs[3])
        elseif reac == 2 # A unravels
            return(rs[1],rs[1]+rs[4])
        elseif reac == 3 # A decays
            return(rs[4],rs[1]+rs[4])
        elseif reac == 4 # A regenerated
            return(rs[3],rs[2]+rs[3])
        elseif reac == 5 # B produced
            return(rs[6],rs[6]+rs[7])
        elseif reac == 6 # B unravels
            return(rs[5],rs[5]+rs[8])
        elseif reac == 7 # B decays
            return(rs[7],rs[5]+rs[8])
        else # B regenerated
            return(rs[8],rs[6]+rs[7])
        end
    else
        error("Invalid reaction code returned")
    end
end

# function to now run a gillespie simulation of the two states
function gillespie(steadA::Array{Float64,1},steadB::Array{Float64,1},mid::Array{Float64,1},ps::Array{Float64,1},
                    noits::Int64,Ω::Int64,pf::Array{Float64,1},pft::Array{Float64,1},pb::Array{Float64,1},
                    pbt::Array{Float64,1})
    # extract all parameters and then rescale
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    k = ps[1]*Ω
    kmin = ps[2]
    q = ps[3]*Ω
    qmin = ps[4]
    K = ps[5]
    Kmin = ps[6]*Ω
    Q = ps[7]
    Qmin = ps[8]*Ω
    r = ps[9]
    f = ps[10]/(Ω^2)
    # set up these gloablly
    times = zeros(2)
    vars = fill(0,2,2)
    SA = SB = SAt = SBt = 0
    # setup for finding if left steady state
    sad = round.(mid*Ω)
    fin = false
    while fin == false
        crossA = false
        crossB = false
        times[1] = 0
        vars[:,1] = round.(steadA*Ω)
        reac = 0
        minA = vars[1,1]
        maxA = vars[1,1]
        minB = vars[2,1]
        maxB = vars[2,1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            if i != 1
                # now use reac to calculate reverse rate
                pb[i-1], pbt[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[:,2], pf[i], pft[i], reac = step(rs,vars[:,1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
                pb[end], pbt[end] = rev(rs,reac)
            end
            # check if max and minimums should be updated
            if vars[1,2] < minA
                minA = vars[1,2]
            elseif vars[1,2] > maxA
                maxA = vars[1,2]
            elseif vars[2,2] < minB
                minB = vars[2,2]
            elseif vars[2,2] > maxB
                maxB = vars[2,2]
            end
            vars[:,1] = vars[:,2]
            times[1] = times[2]
        end
        SA = 0
        for i = 1:length(pb)
            SA += log(pf[i]) - log(pb[i])
        end
        SA = SA/Ω
        # Repeat process for the apparent entropy production
        SAt = 0
        for i = 1:length(pbt)
            SAt += log(pft[i]) - log(pbt[i])
        end
        SAt = SAt/Ω
        # check simulated trajectories remain in bounds
        if (steadA[1] - mid[1]) < 0
            if sad[1] < maxA
                crossA = true
            end
        else
            if sad[1] > minA
                crossA = true
            end
        end
        if (steadA[2] - mid[2]) < 0
            if sad[2] < maxB
                crossB = true
            end
        else
            if sad[2] > minB
                crossB = true
            end
        end
        if crossA == false || crossB == false
            fin = true
            break
        else
            println(steadA)
            println(mid)
            flush(stdout)
        end
    end

    # Second run for other state
    fin = false
    while fin == false
        crossA = false
        crossB = false
        times[1] = 0
        vars[:,1] = round.(steadB*Ω)
        reac = 0
        minA = vars[1,1]
        maxA = vars[1,1]
        minB = vars[2,1]
        maxB = vars[2,1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            if i != 1
                # now use reac to calculate reverse rate
                pb[i-1], pbt[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[:,2], pf[i], pft[i], reac = step(rs,vars[:,1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
                pb[end], pbt[end] = rev(rs,reac)
            end
            # check if max and minimums should be updated
            if vars[1,2] < minA
                minA = vars[1,2]
            elseif vars[1,2] > maxA
                maxA = vars[1,2]
            elseif vars[2,2] < minB
                minB = vars[2,2]
            elseif vars[2,2] > maxB
                maxB = vars[2,2]
            end
            vars[:,1] = vars[:,2]
            times[1] = times[2]
        end
        SB = 0
        for i = 1:length(pb)
            SB += log(pf[i]) - log(pb[i])
        end
        SB = SB/Ω
        # Repeat for reduced master equation case
        SBt = 0
        for i = 1:length(pbt)
            SBt += log(pft[i]) - log(pbt[i])
        end
        SBt = SBt/Ω
        # check simulated trajectories remain in bounds
        if (steadB[1] - mid[1]) < 0
            if sad[1] < maxA
                crossA = true
            end
        else
            if sad[1] > minA
                crossA = true
            end
        end
        if (steadB[2] - mid[2]) < 0
            if sad[2] < maxB
                crossB = true
            end
        else
            if sad[2] > minB
                crossB = true
            end
        end
        if crossA == false || crossB == false
            fin = true
            break
        else
            println(steadB)
            println(mid)
            flush(stdout)
        end
    end
    return(SA,SAt,SB,SBt)
end

# function to now run a gillespie simulation of the two states
function gillespie(steadA::Array{Float64,1},steadB::Array{Float64,1},mid::Array{Float64,1},ps::Array{Float64,1},
                    noits::Int64,Ω::Int64,trajA::Array{Int64,2},trajB::Array{Int64,2},tA::Array{Float64,1},tB::Array{Float64,1})
    # extract all parameters and then rescale
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    k = ps[1]*Ω
    kmin = ps[2]
    q = ps[3]*Ω
    qmin = ps[4]
    K = ps[5]
    Kmin = ps[6]*Ω
    Q = ps[7]
    Qmin = ps[8]*Ω
    r = ps[9]
    f = ps[10]/(Ω^2)
    # set up these gloablly
    times = zeros(2)
    vars = fill(0,2,2)
    # setup for finding if left steady state
    sad = round.(mid*Ω)
    fin = false
    while fin == false
        crossA = false
        crossB = false
        times[1] = 0
        vars[:,1] = round.(steadA*Ω)
        trajA[1,:] = vars[:,1]
        minA = vars[1,1]
        maxA = vars[1,1]
        minB = vars[2,1]
        maxB = vars[2,1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            tA[i+1] = times[2]
            # do gillepsie step
            vars[:,2] = step(rs,vars[:,1])
            trajA[i+1,:] = vars[:,2]
            # check if max and minimums should be updated
            if vars[1,2] < minA
                minA = vars[1,2]
            elseif vars[1,2] > maxA
                maxA = vars[1,2]
            elseif vars[2,2] < minB
                minB = vars[2,2]
            elseif vars[2,2] > maxB
                maxB = vars[2,2]
            end
            vars[:,1] = vars[:,2]
            times[1] = times[2]
        end
        # check simulated trajectories remain in bounds
        if (steadA[1] - mid[1]) < 0
            if sad[1] < maxA
                crossA = true
            end
        else
            if sad[1] > minA
                crossA = true
            end
        end
        if (steadA[2] - mid[2]) < 0
            if sad[2] < maxB
                crossB = true
            end
        else
            if sad[2] > minB
                crossB = true
            end
        end
        if crossA == false || crossB == false
            fin = true
            break
        else
            println(steadA)
            println(mid)
            flush(stdout)
        end
    end

    # Second run for other state
    fin = false
    while fin == false
        crossA = false
        crossB = false
        times[1] = 0
        vars[:,1] = round.(steadB*Ω)
        trajB[1,:] = vars[:,1]
        minA = vars[1,1]
        maxA = vars[1,1]
        minB = vars[2,1]
        maxB = vars[2,1]
        for i = 1:noits
            # calculate rates
            rs = rates(vars[1,1],vars[2,1],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            tB[i+1] = times[2]
            # do gillepsie step
            vars[:,2] = step(rs,vars[:,1])
            trajB[i+1,:] = vars[:,2]
            # check if max and minimums should be updated
            if vars[1,2] < minA
                minA = vars[1,2]
            elseif vars[1,2] > maxA
                maxA = vars[1,2]
            elseif vars[2,2] < minB
                minB = vars[2,2]
            elseif vars[2,2] > maxB
                maxB = vars[2,2]
            end
            vars[:,1] = vars[:,2]
            times[1] = times[2]
        end
        # check simulated trajectories remain in bounds
        if (steadB[1] - mid[1]) < 0
            if sad[1] < maxA
                crossA = true
            end
        else
            if sad[1] > minA
                crossA = true
            end
        end
        if (steadB[2] - mid[2]) < 0
            if sad[2] < maxB
                crossB = true
            end
        else
            if sad[2] > minB
                crossB = true
            end
        end
        if crossA == false || crossB == false
            fin = true
            break
        else
            println(steadB)
            println(mid)
            flush(stdout)
        end
    end
    return(trajA,trajB,tA,tB)
end

# Main function
function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/SupData/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    len = countlines(infile)
    w = 10
    ps = zeros(len,w)
    open(infile, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,w+1)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            for i = 1:w
                ps[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Check there is a file of steady states to be read
    infile = "../Results/SupData/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    len = countlines(infile)
    w = 6
    stead = zeros(len,w)
    open(infile, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,w+1)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            for i = 1:w
                stead[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    # Now do gillespie simulations to find entropy productions
    N = 5000 # number trajectories
    Ω = 5000
    f = [0.0, 0.0]
    noits = 100000
    SLangA = zeros(N,len)
    SLangB = zeros(N,len)
    SMastA = zeros(N,len)
    SMastB = zeros(N,len)
    SRedMA = zeros(N,len)
    SRedMB = zeros(N,len)
    pf = zeros(noits)
    pb = zeros(noits)
    pft = zeros(noits)
    pbt = zeros(noits)
    trajA = zeros(Int64,noits+1,2)
    trajB = zeros(Int64,noits+1,2)
    tA = zeros(noits+1)
    tB = zeros(noits+1)
    for i = 1:len
        println("Starting analysis of parameter set $(i)")
        for j = 1:N
            # Keep user informed of progress
            if j % 100 == 0
                println("Simulation $(j) being performed")
            end
            SMastA[j,i], SRedMA[j,i], SMastB[j,i], SRedMB[j,i] = gillespie([stead[i,1],stead[i,2]],
            [stead[i,5],stead[i,6]],[stead[i,3],stead[i,4]],ps[i,:],noits,Ω,pf,pft,pb,pbt)
            trajA, trajB, tA, tB = gillespie([stead[i,1],stead[i,2]],[stead[i,5],stead[i,6]],[stead[i,3],stead[i,4]],ps[i,:],noits,Ω,trajA,trajB,tA,tB)
            SLangA[j,i] = LangEnt(trajA/Ω,ps[i,:],tA)
            SLangB[j,i] = LangEnt(trajB/Ω,ps[i,:],tB)
        end
    end
    # Now write out data on the steady states
    for i = 1:len
        outfile = "../Results/SupData/Stead$(i)$(ARGS[1]).csv"
        out_file = open(outfile, "w")
        for j = 1:N
            line = "$(SLangA[j,i]),$(SLangB[j,i]),$(SMastA[j,i]),$(SMastB[j,i]),$(SRedMA[j,i]),$(SRedMB[j,i])\n"
            write(out_file,line)
        end
        close(out_file)
    end
    return(nothing)
end

function plotting()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/SupData/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now find length of this parameter file
    len = countlines(infile)
    # check that the first input file exists
    infile = "../Results/SupData/Stead1$(ARGS[1]).csv"
    if ~isfile(infile)
        println("Error: 1st relevant data file is missing.")
        return(nothing)
    end
    # use first file to set length of data
    len2 = countlines(infile)
    SLangA = zeros(len,len2)
    SLangB = zeros(len,len2)
    SMastA = zeros(len,len2)
    SMastB = zeros(len,len2)
    SRedMA = zeros(len,len2)
    SRedMB = zeros(len,len2)
    # Make latex strings to be used in plots
    LM = L"\Delta S_{M}"
    LL = L"\Delta S_{L}"
    for i = 1:len
        # should do some checks on data here
        infile = "../Results/SupData/Stead$(i)$(ARGS[1]).csv"
        if ~isfile(infile)
            println("Error: $(i)th relevant data file is missing.")
            return(nothing)
        end
        # use first file to set length of data
        len3 = countlines(infile)
        if len3 != len2
            println("Error: Length of $(i)th relevant data file does not match with 1st.")
            return(nothing)
        end
        # Now read data from this file
        open(infile, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 1
            for line in eachline(in_file)
                # parse line by finding commas
                L = length(line)
                comma = fill(0,7)
                j = 1
                for l = 1:L
                    if line[l] == ','
                        j += 1
                        comma[j] = l
                    end
                end
                comma[end] = L+1
                SLangA[i,k] = parse(Float64,line[(comma[1]+1):(comma[2]-1)])
                SLangB[i,k] = parse(Float64,line[(comma[2]+1):(comma[3]-1)])
                SMastA[i,k] = parse(Float64,line[(comma[3]+1):(comma[4]-1)])
                SMastB[i,k] = parse(Float64,line[(comma[4]+1):(comma[5]-1)])
                SRedMA[i,k] = parse(Float64,line[(comma[5]+1):(comma[6]-1)])
                SRedMB[i,k] = parse(Float64,line[(comma[6]+1):(comma[7]-1)])
                k += 1
            end
        end
        # now stage to find values to make bins
        maxA = minA = maxB = minB = 0
        maxMA = maximum(SMastA)
        maxLA = maximum(SLangA)
        maxRA = maximum(SRedMA)
        if maxLA > maxRA
            maxA = maxLA
        else
            maxA = maxRA
        end
        minMA = minimum(SMastA)
        minLA = minimum(SLangA)
        minRA = maximum(SRedMA)
        if minLA < minRA
            minA = minLA
        else
            minA = minRA
        end
        maxMB = maximum(SMastB)
        maxLB = maximum(SLangB)
        maxRB = maximum(SRedMB)
        if maxLB > maxRB
            maxB = maxLB
        else
            maxB = maxRB
        end
        minMB = minimum(SMastB)
        minLB = minimum(SLangB)
        minRB = minimum(SRedMB)
        if minLB < minRB
            minB = minLB
        else
            minB = minRB
        end
        # Now use to setup bins
        Nbins = 500
        binsA = range(minA,stop=maxA,length=Nbins)
        binsB = range(minB,stop=maxB,length=Nbins)
        binsMA = range(minMA,stop=maxMA,length=Nbins)
        binsMB = range(minMB,stop=maxMB,length=Nbins)
        # Now carry out plotting
        pyplot() # activate pyplot
        histogram(SRedMA[i,:],bins=binsA,label="",title="Reduced Master Eq")
        plot!(xlabel=LM,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)1RedMStead.png")
        histogram(SRedMB[i,:],bins=binsB,label="",title="Reduced Master Eq")
        plot!(xlabel=LM,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)2RedMStead.png")
        histogram(SMastA[i,:],bins=binsMA,label="",title="Underlying Master Eq")
        plot!(xlabel=LM,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)1MastStead.png")
        histogram(SMastB[i,:],bins=binsMB,label="",title="Underlying Master Eq")
        plot!(xlabel=LM,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)2MastStead.png")
        histogram(SLangA[i,:],bins=binsA,label="",title="Langevin")
        plot!(xlabel=LL,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)1LangStead.png")
        histogram(SLangB[i,:],bins=binsB,label="",title="Langevin")
        plot!(xlabel=LL,ylabel="Number of Trajectories",dpi=300,titlefontsize=20,guidefontsize=16,legendfontsize=12)
        savefig("../Results/SupGraphs/$(i)2LangStead.png")
    end
    return(nothing)
end

@time plotting()
# @time main()
