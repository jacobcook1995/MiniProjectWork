#!/usr/bin/env julia
# switchentS.jl
# A script to calculate the entropy production of the switches between states from the
# gillespie simulation of the chemical master equation, this time for the Schlögl model
#
# Author: Jacob Cook
# Date: March 2019

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64)
    rates = [ k1, K1*X, k2*X*(X-1)*(X-2), K2*X*(X-1) ]
    return(rates)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Int64)
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars += 1 # X produced
        reac = 1
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels
        reac = 2
    elseif r < sum(rs[1:3])
        vars -= 1 # X decays
        reac = 3
    else
        vars += 1 # X regenerated
        reac = 4
    end
    return(vars,reac)
end


function gillespie(ps::Array{Float64,1},star::Int64,fin::Int64,Ω::Float64)
    # extract all parameters and then rescale
    # ps = [ k1, K1, k2, K2 ]
    k1 = ps[1]*Ω
    K1 = ps[2]
    k2 = ps[3]/(Ω^2)
    K2 = ps[4]/Ω
    if star > fin
        high = true
    else
        high = false
    end
    # set up arrays to store paths
    L = 100000000
    ABf = zeros(Int64,L)
    reacs = zeros(Int64,L)
    # run until a switching is observed
    swit = false
    ABf0 = copy(star)
    ABf[1] = ABf0
    j = forj = 1
    # t1 = time_ns()
    while swit == false
        j += 1
        # calculate rates
        rs = rates(ABf0,k1,K1,k2,K2)
        # do gillepsie step
        ABf0, reac0 = step(rs,ABf0)
        # add to vectors
        ABf[j] = ABf0
        reacs[j-1] = reac0
        # stopping condition
        if high == true && ABf0 <= fin
            swit = true
            forj = j
        elseif high == false && ABf0 >= fin
            swit = true
            forj = j
        elseif j >= L # overflow condition
            j = 1
            ABf[1] = ABf[end]
            # No need to update reacs here as none remain valid
        end
    end
    # t2 = time_ns()
    # println("Gillespie Done!")
    # println("Time Elapsed: $((t2-t1)*10.0^-9)")
    # flush(stdout)
    # now remove invalid data from array
    ABf2 = ABf[1:forj]
    reacs2 = reacs[1:forj-1]
    fori = 1
    # now search for sensible starting point for trajectories
    for i = forj:(-1):2
        if high == true && ABf2[i] >= star
            fori = i
            break
        elseif high == false && ABf2[i] <= star
            fori = i
            break
        end
    end
    # remove uneccesary earlier trajectory
    ABf3 = ABf2[fori:end]
    reacs3 = reacs2[fori:end] # change maybe
    return(ABf3,reacs3)
end

# function to find probability of forward jump
# [ k1, K1, k2, K2 ]
function prob(reac::Int64,p1::Int64,p2::Int64,ps::Array{Float64,1})
    rs = [ ps[1], ps[2]*p1, ps[3]*p1*(p1-1)*(p1-2), ps[4]*p1*(p1-1) ]
    rs = rs/sum(rs) # rescale rates to probability
    return(rs[reac]) # then return appropriate probability
end

# function to find probability of reverse jump
function rev(reac::Int64,p1::Int64,p2::Int64,ps::Array{Float64,1})
    rs = [ ps[1], ps[2]*p1, ps[3]*p1*(p1-1)*(p1-2), ps[4]*p1*(p1-1) ]
    rs = rs/sum(rs) # rescale rates to probability
    if reac % 2 == 0
        return(rs[reac-1])
    else
        return(rs[reac+1])
    end
end


# function to find the entropy production of a path
function entp(path::Array{Int64,1},reacs::Array{Int64,1},ps::Array{Float64,1},Ω::Float64)
    bpath = path[end:-1:1]
    breacs = reacs[end:-1:1]
    # rescale constants appropriately
    ps[1] = ps[1]*Ω
    ps[3] = ps[3]/(Ω^2)
    ps[4] = ps[4]/Ω
    pf = zeros(length(reacs))
    pb = zeros(length(reacs))
    for i = 1:length(reacs)
        pf[i] = prob(reacs[i],path[i],path[i+1],ps)
        pb[i] = rev(breacs[i],bpath[i],bpath[i+1],ps)
    end
    ent = log.(pf./pb)
    return(ent)
end

# function to find volume where simulation is relatively easy
function selectvol(A::Array{Float64,1})
    Ω = 0.0
    if length(A) != 2
        println("ERROR: should have provided a two elements for the action!")
        return(nothing)
    end
    maxA = maximum(A)
    if maxA < 0.05
        Ω = 7.0/maxA
    elseif maxA > 2.5
        Ω = 11.0/maxA
    else
        Ω = 10.0/maxA
    end
    return(Ω)
end

function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig3DataS/$(ARGS[1])paraS.csv"
    if ~isfile(infile)
        println("Error: No file of parameters to be read.")
        return(nothing)
    end
    # now read in parameters
    l = countlines(infile)
    w = 4
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
    infile = "../Results/Fig3DataS/$(ARGS[1])steadS.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    l = countlines(infile)
    w = 3
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
    # read in actions now
    acts = zeros(l,8)
    datayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2lD.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2hD.csv"
            end
            # check it exits
            if isfile(infile)
                w = 4
                open(infile, "r") do in_file
                    # only one line now
                    for line in eachline(in_file)
                        # parse line by finding commas
                        L = length(line)
                        comma = fill(0,w+1)
                        m = 1
                        for k = 1:L
                            if line[k] == ','
                                m += 1
                                comma[m] = k
                            end
                        end
                        comma[end] = L+1
                        for k = 1:w
                            acts[i,(j-1)*4+k] = parse(Float64,line[(comma[k]+1):(comma[k+1]-1)])
                        end
                    end
                end
            else # if file doesn't exist inform user of missing data and exclude from plots
                if j == 1
                    println("Missing data for run $(i) high A to high B")
                    datayn[i] = false
                else
                    println("Missing data for run $(i) high B to high A")
                    datayn[i] = false
                end
            end
        end
    end
    # Setting parameter sets to look at
    # I = collect(1:100)
    I = [5,10,20,21,36,71,81,99,100] # REMOVE WHEN DONE
    for i = 1:length(I)
        println("Run $(I[i]) Started.")
        flush(stdout)
        # now calulate entropy productions of switching
        # calculate optimum value each time
        Ωs = selectvol([acts[I[i],2],acts[I[i],6]])
        # Need to increase Ωs in the hope it improves result
        Ωs = 1.25*Ωs
        len = 100 # run 100 for each, forward and backward
        entpf = zeros(len)
        entpb = zeros(len)
        star = round.(Int64,Ωs*steads[I[i],1])
        fin = round.(Int64,Ωs*steads[I[i],3])
        for j = 1:len
            path, reacs = gillespie(ps[I[i],:],star,fin,Ωs)
            entprod = entp(path,reacs,ps[I[i],:],Ωs)
            entpf[j] = sum(entprod)/Ωs
            # opposite direction
            path, reacs = gillespie(ps[I[i],:],fin,star,Ωs)
            entprod = entp(path,reacs,ps[I[i],:],Ωs)
            entpb[j] = sum(entprod)/Ωs
            if j % 10 == 0
                println("Step $(j) done.")
                flush(stdout)
            end
        end
        # Then output all data
        # Do A to B first
        outfile = "../Results/Fig3DataS/Traj/$(I[i])$(ARGS[1])h2lMast.csv"
        out_file = open(outfile, "w")
        line = ""
        for j = 1:len
            line *= "$(entpf[j]),"
        end
        # add a new line and volume
        line *= "$(Ωs)\n"
        write(out_file,line)
        close(out_file)
        # The do B to A
        outfile = "../Results/Fig3DataS/Traj/$(I[i])$(ARGS[1])l2hMast.csv"
        out_file = open(outfile, "w")
        line = ""
        for j = 1:len
            line *= "$(entpb[j]),"
        end
        # add a new line and volume
        line *= "$(Ωs)\n"
        write(out_file,line)
        close(out_file)
    end
    return(nothing)
end

@time main()
