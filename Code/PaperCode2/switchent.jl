#!/usr/bin/env julia
# switchent.jl
# A script to calculate the entropy production of the switches between states from the
# gillespie simulation of the chemical master equation
#
# Author: Jacob Cook
# Date: December 2018

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to advance gillespie one step
function step(rates::Array{Float64,1},vars::Array{Int64,1})
    r = rand()
    rs = rates/sum(rates)
    if r < rs[1]
        vars[1] += 1 # A produced
        reac = 1
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
        reac = 3
    elseif r < sum(rs[1:4])
        vars[1] += 1 # A regenerated
        reac = 4
    elseif r < sum(rs[1:5])
        vars[2] += 1 # B produced
        reac = 5
    elseif r < sum(rs[1:6])
        vars[2] -= 1 # B unravels
        reac = 6
    elseif r < sum(rs[1:7])
        vars[2] -= 1 # B decays
        reac = 7
    else
        vars[2] += 1 # B regenerated
        reac = 8
    end
    return(vars,reac)
end

# function to actually run a gillepsie simulation
function gillespie(ps::Array{Float64,1},star::Array{Int64,1},fin::Array{Int64,1},Ω::Float64)
    k = ps[1]
    kmin = ps[2]
    q = ps[3]
    qmin = ps[4]
    K = ps[5]
    Kmin = ps[6]
    Q = ps[7]
    Qmin = ps[8]
    r = ps[9]
    f = ps[10]
    if star[1] > fin[1]
        highA = true
    else
        highA = false
    end
    # rescale constants appropraitly
    Kmin = Kmin*Ω
    Qmin = Qmin*Ω
    k = k*Ω
    q = q*Ω
    f = f/(Ω^2)
    # set up arrarys to store paths
    L = 100000000
    ABf = zeros(Int64,2,L)
    reacs = zeros(Int64,L)
    # run until a switching is observed
    swit = false
    ABf0 = copy(star)
    ABf[:,1] = ABf0
    j = forj = 1
    # t1 = time_ns()
    while swit == false
        j += 1
        # calculate rates
        rs = rates(ABf0[1],ABf0[2],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        # do gillepsie step
        ABf0, reac0 = step(rs,ABf0)
        # add to vectors
        ABf[:,j] = ABf0
        reacs[j-1] = reac0
        # stopping condition
        if highA == true && ABf0[1] <= fin[1] && ABf0[2] >= fin[2]
            swit = true
            forj = j
        elseif highA == false && ABf0[1] >= fin[1] && ABf0[2] <= fin[2]
            swit = true
            forj = j
        elseif j >= L # overflow condition
            j = 1
            ABf[:,1] = ABf[:,end]
            # No need to update reacs here as none remain valid
        end
    end
    # t2 = time_ns()
    # println("Gillespie Done!")
    # println("Time Elapsed: $((t2-t1)*10.0^-9)")
    # flush(stdout)
    # now remove invalid data from array
    ABf2 = ABf[:,1:forj]
    reacs2 = reacs[1:forj-1]
    fori = 1
    # now search for sensible starting point for trajectories
    for i = forj:(-1):2
        if highA == true && ABf2[1,i] >= star[1] && ABf2[2,i] <= star[2]
            fori = i
            break
        elseif highA == false && ABf2[1,i] <= star[1] && ABf2[2,i] >= star[2]
            fori = i
            break
        end
    end
    # remove uneccesary earlier trajectory
    ABf3 = ABf2[:,fori:end]
    reacs3 = reacs2[fori:end] # change maybe
    return(ABf3,reacs3)
end

# function to find probability of forward jump
# [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
function prob(reac::Int64,p1::Array{Int64,1},p2::Array{Int64,1},ps::Array{Float64,1})
    rs = [ ps[9]*ps[1]/(ps[9]+ps[10]*p1[2]*(p1[2]-1)), ps[2]*p1[1], ps[5]*p1[1], ps[6], ps[9]*ps[3]/(ps[9]+ps[10]*p1[1]*(p1[1]-1)),
    ps[4]*p1[2], ps[7]*p1[2], ps[8] ]
    rs = rs/sum(rs) # rescale rates to probability
    return(rs[reac]) # then return appropriate probability
end

# function to find probability of reverse jump
function rev(reac::Int64,p1::Array{Int64,1},p2::Array{Int64,1},ps::Array{Float64,1})
    rs = [ ps[9]*ps[1]/(ps[9]+ps[10]*p1[2]*(p1[2]-1)), ps[2]*p1[1], ps[5]*p1[1], ps[6], ps[9]*ps[3]/(ps[9]+ps[10]*p1[1]*(p1[1]-1)),
    ps[4]*p1[2], ps[7]*p1[2], ps[8] ]
    rs = rs/sum(rs) # rescale rates to probability
    if reac % 2 == 0
        return(rs[reac-1])
    else
        return(rs[reac+1])
    end
end

# function to find the entropy production of a path
function entp(path::Array{Int64,2},reacs::Array{Int64,1},ps::Array{Float64,1},Ω::Float64)
    bpath = path[:,end:-1:1]
    breacs = reacs[end:-1:1]
    # rescale constants appropraitly
    ps[6] = ps[6]*Ω
    ps[8] = ps[8]*Ω
    ps[1] = ps[1]*Ω
    ps[3] = ps[3]*Ω
    ps[10] = ps[10]/(Ω^2)
    pf = zeros(length(reacs))
    pb = zeros(length(reacs))
    for i = 1:length(reacs)
        pf[i] = prob(reacs[i],path[:,i],path[:,i+1],ps)
        pb[i] = rev(breacs[i],bpath[:,i],bpath[:,i+1],ps)
    end
    ent = log.(pf./pb)
    return(ent)
end

function selectvol(A::Array{Float64,1},steads::Array{Float64,1})
    Ω = 0.0
    if length(A) != 2
        println("ERROR: should have provided a two elements for the action!")
        return(nothing)
    end
    if length(steads) != 6
        println("ERROR: should have provided a six points to describe states!")
        return(nothing)
    end
    # Check for differnces in steady state height
    δtot = (steads[1] - steads[6])/(steads[1] + steads[6])
    maxA = maximum(A)
    hA = (steads[1]+steads[2]) # find heights of states
    hB = (steads[5]+steads[6])
    hs = (steads[3]+steads[4])
    if maxA > 20.0 # catch really high values
        Ω = 3.0/maxA
    elseif abs(δtot) > 0.40 # not super high but skewed
        Ω = 2.5/maxA
    elseif maxA > 4.5 || hs/hA < 0.125 || hs/hB < 0.125
        Ω = 4.5/maxA
    elseif hs/hA < 0.25 || hs/hB < 0.25 # low stable state condidtion
        Ω = 5.5/maxA
    else
        Ω = 7.5/maxA
    end
    # step to undo volume reduction when actions are very asymmetric
    δA = abs(A[1]-A[2])
    if δA/maxA >= 0.9 && abs(δtot) < 0.5 && maxA < 8.0
        Ω = 7.5/maxA
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
    acts = zeros(l,8)
    datayn = fill(true,l)
    for i = 1:l
        for j = 1:2
            # set file to read in
            if j == 1
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])A2BD.csv"
            else
                infile = "../Results/Fig3Data/Traj/$(i)$(ARGS[1])B2AD.csv"
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
    # Full parameter set
    I = vcat(collect(5:10),12,collect(14:18),collect(20:31),collect(33:47),collect(49:100))
    for i = 1:length(I)
        println("Run $(I[i]) Started.")
        flush(stdout)
        # now calulate entropy productions of switching
        # calculate optimum value each time
        Ωs = selectvol([acts[I[i],2],acts[I[i],6]],steads[I[i],:])
        len = 100 # run 100 for each, forward and backward
        entpf = zeros(len)
        entpb = zeros(len)
        star = round.(Int64,Ωs*[steads[I[i],1],steads[I[i],2]])
        fin = round.(Int64,Ωs*[steads[I[i],5],steads[I[i],6]])
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
        outfile = "../Results/Fig3Data/Traj/$(I[i])$(ARGS[1])A2BMast.csv"
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
        outfile = "../Results/Fig3Data/Traj/$(I[i])$(ARGS[1])B2AMast.csv"
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
