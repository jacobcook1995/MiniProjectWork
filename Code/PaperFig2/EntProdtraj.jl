#!/usr/bin/env julia
# EntProdtraj.jl
# A script to read in parameters and trajcetories and then to generate nice plots
# of the entropy productions by two differenet methods.
# 1) using the langevin entropy production term
# 2) by using a large volume reduced form of the master equation
#
# Author: Jacob Cook
# Date: January 2019

using SymEngine
using Plots
using LaTeXStrings
using PyCall
import PyPlot

# make a symbolic diffusion matrix
function Ds()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,2,2)
    e[1,2] = e[2,1] = 0
    e[1,1] = sqrt(k*r/(r + f*B^2) + kmin*A + K*A + Kmin)
    e[2,2] = sqrt(q*r/(r + f*A^2) + qmin*B + Q*B + Qmin)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    for j = 1:2
        for i = 1:2
            Dmin[i,j] = expand(Dmin[i,j])
        end
    end
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,2)
    b[1] = k*r/(r + f*B^2) - kmin*A - K*A + Kmin
    b[2] = q*r/(r + f*A^2) - qmin*B - Q*B + Qmin
    return(b)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y1, y2 = symbols("y1 y2")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1] + b[1]*Dmin[1,2]*b[2]
    num += b[2]*Dmin[2,1]*b[1] + b[2]*Dmin[2,2]*b[2]
    den = 0
    den += y1*Dmin[1,1]*y1 + y1*Dmin[1,2]*y2
    den += y2*Dmin[2,1]*y1 + y2*Dmin[2,2]*y2
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2 = symbols("y1 y2")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,2)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    ϑ = Array{SymEngine.Basic,1}(undef,2)
    ϑ[1] = Dmin[1,1]*c[1] + Dmin[1,2]*c[2]
    ϑ[2] = Dmin[2,1]*c[1] + Dmin[2,2]*c[2]
    return(ϑ)
end

# function to generate one of each symbolic object and to substitute parameter values into it
# this is a computationally intensive function, need to minimize the number of calls
function gensyms(ps::Array{Float64,1})
    # create symbolic objects
    ϑ = ϑs()
    λ = λs()
    b = bs()
    Dmin = Dmins()
    # specify symbols that will be substituted for
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("K k Q q kmin Kmin qmin Qmin f r")
    # now perform substitutions
    λ = subs(λ, k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
    λ = subs(λ, Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    for i = 1:2
        ϑ[i] = subs(ϑ[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        ϑ[i] = subs(ϑ[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        b[i] = subs(b[i], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
        b[i] = subs(b[i], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
    end
    for j = 1:2
        for i = 1:2
            Dmin[i,j] = subs(Dmin[i,j], k=>ps[1], kmin=>ps[2], q=>ps[3], qmin=>ps[4], K=>ps[5], Kmin=>ps[6])
            Dmin[i,j] = subs(Dmin[i,j], Q=>ps[7], Qmin=>ps[8], r=>ps[9], f=>ps[10])
        end
    end
    return(ϑ,λ,b,Dmin)
end

# function to generate the variables needed for a given path
function genvars(x::Array{Float64,2},λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int64,Nmid::Int64)
    # define neccesary symbols
    A, B, y1, y2 = symbols("A B y1 y2")
    # calculate velocities
    xprim = fill(NaN, NG+1, 2)
    for i = 2:NG
        for j = 1:2
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 2)
    ϑt = Array{SymEngine.Basic,1}(undef,2)
    for j = 1:2
        for i = 2:NG
            ϑt[j] = subs(ϑ[j], A=>x[i,1], B=>x[i,2])
            ϑs[i,j] = subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2]) |> float
        end
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# function to find the times of each point
function times(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},NG::Int64)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::Array{Float64,1},x::Array{Float64,2},NG::Int64,NM::Int64)
    # Make discrete vector of time points
    t = zeros(NM+1)
    for i = 1:NM+1
        t[i] = (i-1)*ts[end]/NM
    end
    # Find index of first element greater than t[i] in ts
    inds = fill(0,NM+1)
    j = 1
    for i = 1:NM+1
        higher = false
        while higher == false
            if ts[j] >= t[i] || j == NG+1
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they are fixed
    path = zeros(NM+1,2)
    path[1,:] = x[1,:]
    path[NM+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        for j = 1:2
            x₀ = x[one,j]
            x₁ = x[two,j]
            path[i,j] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
        end
    end
    return(path)
end

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
    f[1] = ps[1]*ps[9]/(ps[9]+ps[10]*x[2]*x[2]) - ps[2]*x[1] - ps[5]*x[1] + ps[6]
    f[2] = ps[3]*ps[9]/(ps[9]+ps[10]*x[1]*x[1]) - ps[4]*x[2] - ps[7]*x[2] + ps[8]
    return(f)
end

# A function to find and plot the langevin entropy productions of the various trajectories
function LangEnt(traj::Array{Float64,2},ps::Array{Float64,1},dt::Float64)
    entp = zeros(size(traj,1)-1)
    f = [0.0; 0.0]
    D = [0.0 0.0; 0.0 0.0]
    for i = 1:size(entp,1)
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
    return(entp)
end

# function that takes in two points and finds the probability that the switch that generates them happens
function probs(star::Array{Int64,1},fin::Array{Int64,1},k::Float64,K::Float64,q::Float64,Q::Float64,
                kmin::Float64,Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ k*r/(r+f*star[2]^2), kmin*star[1], K*star[1], Kmin, q*r/(r+f*star[1]^2), qmin*star[2], Q*star[2], Qmin ]
    # not actually possible to distinguish two processes, which will effect the probabilities
    if fin[1] - star[1] == 0
        if fin[2] - star[2] == 1
            P = (rates[5] + rates[8])/sum(rates)
        elseif fin[2] - star[2] == -1
            P = (rates[6] + rates[7])/sum(rates)
        else
            error()
        end
    elseif fin[2] - star[2] == 0
        if fin[1] - star[1] == 1
            P = (rates[1] + rates[4])/sum(rates)
        elseif fin[1] - star[1] == -1
            P = (rates[2] + rates[3])/sum(rates)
        else
            error()
        end
    else
        println("Error: Appears step is in both directions")
        error()
    end
    return(P)
end

# A function to find and plot the reduced master equation entropy productions of the various trajectories
function MastEnt(traj::Array{Float64,2},ps::Array{Float64,1},T::Float64,Ω::Int64)
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
    # rescale rates appropriately
    k = k*Ω
    Kmin = Kmin*Ω
    q = q*Ω
    Qmin = Qmin*Ω
    f = f/(Ω^2)
    ps = [k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    # rescale trajectory
    traj = traj*Ω
    # Should now make path in A,B
    path = Array{Int64,2}(undef,0,2)
    ts = Array{Float64,1}(undef,0)
    # approximate steady state
    val = fill(0,1,2)
    for i = 1:2
        val[i] = round(Int64,traj[1,i])
    end
    # initial time t=0
    t = 0
    ts = vcat(ts,t)
    # and add as start of path
    path = vcat(path,val)
    L = size(traj,1)
    # want to find point where A/B crossed the 0.5 point
    for i = 2:L
        # first establish if there has been any change from the prior step
        if round(Int64,traj[i,1]) == path[end,1] && round(Int64,traj[i,2]) == path[end,2]
            # if no change between successive points no reason to update path
        else
            # count changes to A and B between two points
            dA = round(Int64,traj[i,1]) - path[end,1]
            dB = round(Int64,traj[i,2]) - path[end,2]
            th = (i-1)*T/(L-1) # Time here
            tp = (i-2)*T/(L-1) # Time at prior step
            vals = Array{Int64,2}(undef,0,2) # array to store values
            tempt = Array{Float64,1}(undef,0) # vector to temporaily store times
            if dA != 0
                # find difference in A in interval
                da = (traj[i,1] - traj[i-1,1])
                # find first & potentially only point of change
                ap = round(Int64,traj[i,1]) - 0.5*sign(dA)
                # find distance from point of point of change
                deltaA = traj[i,1] - ap
                # as a fraction of spacing
                fracA = deltaA/da
                # find time and save
                t = th - fracA*(th - tp)
                tempt = vcat(tempt,t)
                # find B and save point
                B = traj[i,2] - fracA*(traj[i,2] - traj[i-1,2])
                val[1] = round(Int64,traj[i,1])
                val[2] = round(Int64,B)
                vals = vcat(vals,val)
                n = 1 # first point done
                # loop until all points are done
                done = false
                while done == false
                    if n == abs(dA)
                        done = true
                        break
                    end
                    deltaA = 1*sign(dA)
                    fracA = deltaA/da
                    # remove time of gap from time
                    tn = fracA*(th - tp)
                    t -= tn
                    tempt = vcat(tempt,t)
                    B -= fracA*(traj[i,2] - traj[i-1,2])
                    val[1] = round(Int64,traj[i,1]) - n*sign(dA)
                    val[2] = round(Int64,B)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            if dB != 0
                # find difference in B in interval
                db = (traj[i,2] - traj[i-1,2])
                # find first & potentially only point of change
                bp = round(Int64,traj[i,2]) - 0.5*sign(dB)
                # find distance from point of point of change
                deltaB = traj[i,2] - bp
                # as a fraction of spacing
                fracB = deltaB/db
                # find time and save
                t = th - fracB*(th - tp)
                tempt = vcat(tempt,t)
                # find A and save point
                A = traj[i,1] - fracB*(traj[i,1] - traj[i-1,1])
                val[1] = round(Int64,A)
                val[2] = round(Int64,traj[i,2])
                vals = vcat(vals,val)
                n = 1 # first point done
                # loop until all points are done
                done = false
                while done == false
                    if n == abs(dB)
                        done = true
                        break
                    end
                    deltaB = 1*sign(dB)
                    fracB = deltaB/db
                    # remove time of gap from time
                    tn = fracB*(th - tp)
                    t -= tn
                    tempt = vcat(tempt,t)
                    A -= fracB*(traj[i,1] - traj[i-1,1])
                    val[1] = round(Int64,A)
                    val[2] = round(Int64,traj[i,2]) - n*sign(dB)
                    vals = vcat(vals,val)
                    n += 1 # another point completed
                end
            end
            # Now reorder vals by the times
            p = sortperm(tempt)
            tempt = sort(tempt)
            vals = vals[p,:]
            # then vcat to path
            ts = vcat(ts,tempt)
            path = vcat(path,vals)
        end
    end
    # need both forward and backward trajectories in order to establish entropy production
    backp = path[end:-1:1,:]
    backts = ts[end:-1:1]
    len = length(ts)-1
    Pf = zeros(len)
    Pb = zeros(len)
    for i = 2:len+1
        Pf[i-1] = probs(path[i-1,:],path[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        Pb[i-1] = probs(backp[i-1,:],backp[i,:],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
    end
    entp = log.(Pf./Pb[end:-1:1])
    # Now should rebin the entropy data into the appropriate segments
    len2 = size(traj,1)-1
    entsr = zeros(len2)
    inds = zeros(Int64,len2)
    # find index to sum up to for each one
    for i = 1:len2-1
        point = round.(Int64,traj[i+1,:])
        ind = findall((path[:,1].==point[1]) .& (path[:,2].==point[2]))
        if length(ind) == 1
            if ind[1] == 1
                inds[i] = 1
            else
                inds[i] = ind[1] - 1
            end
        else
            println("Seems like the path is crossing itself")
            error()
        end
    end
    inds[end] = length(entp)
    # then sum entropies into appropraite sections
    entsr[1] = sum(entp[1:inds[1]])
    for i = 1:len2-1
        entsr[i+1] = sum(entp[inds[i]:inds[i+1]])
    end
    # finally extract the path in A
    entp = entsr
    return(entp)
end

# main function
function main()
    println("Compiled, Starting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Check there is a file of parameters to be read
    infile = "../Results/Fig2Data/$(ARGS[1])para.csv"
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
    # now read in the 6 trajcetories
    N = 600
    traj = zeros(N+1,2*2*len)
    for i = 1:len
        for j = 1:2
            if j == 1
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])A2B.csv"
            else
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])B2A.csv"
            end
            w = 2
            open(infile, "r") do in_file
                # Use a for loop to process the rows in the input file one-by-one
                k = 1
                for line in eachline(in_file)
                    # parse line by finding commas
                    L = length(line)
                    comma = fill(0,w+1)
                    m = 1
                    for l = 1:L
                        if line[l] == ','
                            m += 1
                            comma[m] = l
                        end
                    end
                    comma[end] = L+1
                    for l = 1:w
                        traj[k,4*(i-1)+2*(j-1)+l] = parse(Float64,line[(comma[l]+1):(comma[l+1]-1)])
                    end
                    k += 1
                end
            end
        end
    end
    # These trajectories now must be converted to time discretisation
    NG = NM = N
    traj2 = zeros(NM+1,2*2*len)
    # define NG and Nmid and use to find variables
    Nmid = convert(Int64, ceil((NG+1)/2))
    for i = 1:len
        for j = 1:2
            # So first generate symbolic objects
            ϑ, λ, b, Dmin = gensyms(ps[i,:])
            d = 2*(j-1)+4*(i-1)
            x, xprim, λs, ϑs, λprim = genvars(traj[:,(1+d):(2+d)],λ,ϑ,NG,Nmid)
            # use function Ŝ to find the action associated with this path
            λs[1] = λs[2]
            λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
            λs[end] = λs[end-1]
            # find and save action from geometric method
            tims2 = times(x,xprim,λs,ϑs,λprim,NG)
            # save time of path
            Tp = tims2[end]
            traj2[:,(1+d):(2+d)] = timdis(tims2,x,NG,NM)
        end
    end
    # Now read in additional data from the other stuff
    w = 4
    Act = zeros(2*len,4)
    for i = 1:len
        for j = 1:2
            if j == 1
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])A2BD.csv"
            else
                infile = "../Results/Fig2Data/$(i)$(ARGS[1])B2AD.csv"
            end
            open(infile, "r") do in_file
                # Single line file so doesn't matter
                for line in eachline(in_file)
                    # parse line by finding commas
                    L = length(line)
                    comma = fill(0,w+1)
                    m = 1
                    for l = 1:L
                        if line[l] == ','
                            m += 1
                            comma[m] = l
                        end
                    end
                    comma[end] = L+1
                    for l = 1:w
                        Act[2*(i-1)+j,l] = parse(Float64,line[(comma[l]+1):(comma[l+1]-1)])
                    end
                end
            end
        end
    end
    # Check there is a file of steady states to be read
    infile = "../Results/Fig2Data/$(ARGS[1])stead.csv"
    if ~isfile(infile)
        println("Error: No file of steady states to be read.")
        return(nothing)
    end
    # now read in steady states
    len = countlines(infile)
    w = 6
    steads = zeros(len,w)
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
    # Now calculate and plot langevin entropy productions
    LatS = L"\Delta S_{L}"
    LatS2 = L"\Delta S_{M}"
    pyplot()
    Ω = 5000
    for i = 1:len
        for j = 1:2
            d = 2*(j-1)+4*(i-1)
            entp = LangEnt(traj2[:,(1+d):(2+d)],ps[i,:],Act[2*(i-1)+j,1]/N)
            println("$(i),$(j),$(sum(entp)),$(Act[2*(i-1)+j,4])")
            entp2 = MastEnt(traj2[:,(1+d):(2+d)],ps[i,:],Act[2*(i-1)+j,1],Ω)
            println("$(i),$(j),$(sum(entp2)/Ω)")
            plot(traj2[1:end-1,1+d],entp,label=LatS,dpi=300,legend=:best,title="Langevin")
            plot!(xlabel="Concentration A",ylabel=LatS,titlefontsize=20,guidefontsize=16,legendfontsize=12)
            scatter!([steads[i,1+4*(j-1)]],[0.0],markersize=6,color=:black,label="Start")
            scatter!([steads[i,3]],[0.0],markersize=5,color=:black,markershape=:x,label="Saddle")
            scatter!([steads[i,5-4*(j-1)]],[0.0],markersize=6,color=:white,label="End")
            savefig("../Results/Fig2Graphs/$(i)$(j)LangEnt.png")
            plot(traj2[1:end-1,1+d],entp2/(Ω),label=LatS2,dpi=300,legend=:best,title="Reduced Master Eq")
            plot!(xlabel="Concentration A",ylabel=LatS2,titlefontsize=20,guidefontsize=16,legendfontsize=12)
            scatter!([steads[i,1+4*(j-1)]],[0.0],markersize=6,color=:black,label="Start")
            scatter!([steads[i,3]],[0.0],markersize=5,color=:black,markershape=:x,label="Saddle")
            scatter!([steads[i,5-4*(j-1)]],[0.0],markersize=6,color=:white,label="End")
            savefig("../Results/Fig2Graphs/$(i)$(j)MastEnt.png")
        end
    end
    return(nothing)
end

@time main()
