#!/usr/bin/env julia
# entgill.jl
# A script to read in parameters and steady states and then do Gillespie simulationa
# in order to obtain values of entropy production based on master equation dynamics
#
# Author: Jacob Cook
# Date: November 2018

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
function rev(rs::Array{Float64,1},reac::Int64)
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

function gillespie(steadA::Array{Float64,1},steadB::Array{Float64,1},mid::Array{Float64,1},ps::Array{Float64,1},
                    noits::Int64,Ω::Int64,pf::Array{Float64,1},pb::Array{Float64,1})
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
    SA = SB = 0
    # setup for finding if left steady state
    sad = round.(mid*Ω)
    fin = false
    # println(steadA)
    # println(steadB)
    # println(mid)
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
                pb[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[:,2], pf[i], reac = step(rs,vars[:,1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
                pb[end] = rev(rs,reac)
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
        SA = SA/(Ω*times[2])
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
                pb[i-1] = rev(rs,reac)
            end
            # calculate timestep
            τ = timstep(rs)
            # update time
            times[2] = times[1] + τ
            # do gillepsie step
            vars[:,2], pf[i], reac = step(rs,vars[:,1],reac)
            # final reverse rate
            if i == noits
                rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
                pb[end] = rev(rs,reac)
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
        SB = SB/(Ω*times[2])
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
    return([SA,SB])
end

function main()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3Data/$(ARGS[1])para.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 10
    ps = zeros(l,w)
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
    infile = "../Results/Fig3Data/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    l = countlines(infile)
    w = 6
    steads = zeros(l,w)
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
                steads[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
            end
            k += 1
        end
    end
    ents = zeros(l,2)
    # things to fiddle with to get to work
    noits = 100000000
    Ω = 5000
    pf = zeros(noits)
    pb = zeros(noits)
    for i = 1:l
        println(i)
        flush(stdout)
        ents[i,:] = gillespie([steads[i,1],steads[i,2]],[steads[i,5],steads[i,6]],[steads[i,3],steads[i,4]],ps[i,:],noits,Ω,pf,pb)
    end
    # Now write out ents to file
    output_file = "../Results/Fig3Data/$(ARGS[1])gill.csv"
    out_file = open(output_file, "w")
    for i = 1:size(ents,1)
        line = ""
        for j = 1:size(ents,2)
            line *= "$(ents[i,j]),"
        end
        # remove surplus ,
        line = line[1:end-1]
        # add a new line
        line *= "\n"
        write(out_file,line)
    end
    return(nothing)
end


@time main()
