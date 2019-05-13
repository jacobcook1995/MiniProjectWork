# !/usr/bin/env julia
# indent.jl
# A script that reads in the identified best and worst fitting trajectories
# Which then attempts to identify common features
#
# Author: Jacob Cook
# Date: May 2019

using SymEngine
using Plots
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
function timdis(ts::Array{Float64,1},x::Array{Float64,2},NG::Int64,NM::Int64,Tmid::Float64)
    Ns = 1
    signed = false
    # Make discrete vector of time points
    t = zeros(NM+1)
    for i = 1:NM+1
        t[i] = (i-1)*ts[end]/NM
        if signed == false && t[i] >= Tmid
            Ns = i
            signed = true
        end
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
    return(path,Ns)
end

# function to calculate speed along a path
function speed(path::Array{Float64,2},ps::Array{Float64,1})
    v = zeros(size(path,1)-1,2)
    f = zeros(size(path,1)-1,2)
    C = zeros(size(path,1)-1)
    ΔS = zeros(size(path,1)-1)
    Dm = zeros(size(path,1)-1,2,2)
    # Generate values needed for time discretisation
    ϑ, λ, b, Dmin = gensyms(ps)
    # define NG and Nmid and use to find variables
    NG = NM = size(path,1)-1
    Nmid = convert(Int64, ceil((NG+1)/2))
    x, xprim, λs, ϑs, λprim = genvars(path,λ,ϑ,NG,Nmid)
    # fix boundary λ values
    λs[1] = λs[2]
    λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
    λs[end] = λs[end-1]
    # Find times of each point
    tims2 = times(x,xprim,λs,ϑs,λprim,NG)
    # save time of path
    Tp = tims2[end]
    Tmid = tims2[Nmid]
    path2, Ns = timdis(tims2,x,NG,NM,Tmid)
    δt = tims2[2] - tims2[1]
    for i = 1:size(v,1)
        v[i,1] = (path2[i+1,1]-path2[i,1])/δt
        v[i,2] = (path2[i+1,2]-path2[i,2])/δt
    end
    A, B = symbols("A B")
    # Adding calculation of f in
    for i = 1:size(f,1)
        posA = (path2[i+1,1] + path2[i,1])/2
        posB = (path2[i+1,2] + path2[i,2])/2
        bt = copy(b)
        for j = 1:2
            bt[j] = subs(bt[j], A=>posA, B=>posB) |> float
        end
        f[i,1] = bt[1]
        f[i,2] = bt[2]
    end
    # Now calculate diffusion matrices
    for i = 1:size(Dm,1)
        posA = (path2[i+1,1] + path2[i,1])/2
        posB = (path2[i+1,2] + path2[i,2])/2
        Dmt = copy(Dmin)
        for j = 1:2
            for k = 1:2
                Dmt[j,k] = subs(Dmt[j,k], A=>posA, B=>posB) |> float
            end
        end
        Dm[i,1,1] = Dmt[1,1]
        Dm[i,2,2] = Dmt[2,2]
        Dm[i,1,2] = Dm[i,2,1] = 0.0
    end
    # and then calculate conservative actions
    for i = 1:length(C)
        for j = 1:2
            C[i] += 0.5*(f[i,j]*Dm[i,j,j]*f[i,j] + v[i,j]*Dm[i,j,j]*v[i,j])*δt
        end
    end
    # Calculate entropy production
    for i = 1:length(ΔS)
        for j = 1:2
            ΔS[i] += (f[i,j]*Dm[i,j,j]*v[i,j])*δt
        end
    end
    return(v,path2,f,C,ΔS)
end

function main()
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    end
    # Firstly read in all parameters
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
    # Define input files
    infileb = "../Results/Fig3Data/Best.csv"
    infilew = "../Results/Fig3Data/Worst.csv"
    # Check that both exist
    if ~isfile(infileb)
        println("Error: No file of best indices to be read.")
        return(nothing)
    end
    if ~isfile(infilew)
        println("Error: No file of worst indices to be read.")
        return(nothing)
    end
    # Then read both in
    l = countlines(infileb)
    best = zeros(Int64,l)
    open(infileb, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,3)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            # Parse from beginning (comma[1]) to first comma (comma[2])
            temp = parse(Float64,line[(comma[1]+1):(comma[2]-1)])
            # Convert to Int
            best[k] = convert(Int64,temp)
            k += 1
        end
    end
    # Now the worst data
    l = countlines(infilew)
    worst = zeros(Int64,l)
    open(infilew, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,3)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            # Parse from beginning (comma[1]) to first comma (comma[2])
            temp = parse(Float64,line[(comma[1]+1):(comma[2]-1)])
            # Covert to Int
            worst[k] = convert(Int64,temp)
            k += 1
        end
    end
    # Now read in best trajectories and plot
    for i = 1:3#length(best)
        infile1 = "../Results/Fig3Data/Traj/$(best[i])$(ARGS[1])A2B.csv"
        infile2 = "../Results/Fig3Data/Traj/$(best[i])$(ARGS[1])B2A.csv"
        # now should read in 1st path
        l = countlines(infile1)
        w = 2
        path1 = zeros(l,w)
        open(infile1, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 1
            for line in eachline(in_file)
                # parse line by finding commas
                L = length(line)
                comma = fill(0,w+1)
                m = 1
                for i = 1:L
                    if line[i] == ','
                        m += 1
                        comma[m] = i
                    end
                end
                comma[end] = L+1
                for i = 1:w
                    path1[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # now should read in 2nd path
        l = countlines(infile2)
        w = 2
        path2 = zeros(l,w)
        open(infile2, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 1
            for line in eachline(in_file)
                # parse line by finding commas
                L = length(line)
                comma = fill(0,w+1)
                m = 1
                for i = 1:L
                    if line[i] == ','
                        m += 1
                        comma[m] = i
                    end
                end
                comma[end] = L+1
                for i = 1:w
                    path2[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Calculate and then plot velocities
        v1, tpath1, f1, C1, ΔS1 = speed(path1,ps[best[i],:])
        v2, tpath2, f2, C2, ΔS2 = speed(path2,ps[best[i],:])
        plot(path1[:,2],path1[:,1])
        plot!(path2[:,2],path2[:,1])
        savefig("../Results/BestWorst/Best$(i).png")
        plot(tpath1[2:end,2],abs.(v1))
        plot!(tpath2[2:end,2],abs.(v2))
        savefig("../Results/BestWorst/VelocityBest$(i).png")
        plot(tpath1[2:end,2],abs.(f1))
        plot!(tpath2[2:end,2],abs.(f2))
        savefig("../Results/BestWorst/ForceBest$(i).png")
        plot(tpath1[2:end,2],C1)
        plot!(tpath2[2:end,2],C2)
        savefig("../Results/BestWorst/ConservBest$(i).png")
        # Make cumaltive actions
        CC1 = zeros(length(C1))
        CC2 = zeros(length(C2))
        for i = 1:length(C1)
            CC1[i] = sum(C1[1:i])
            CC2[i] = sum(C2[end+1-i:end])
        end
        plot(tpath1[2:end,2],CC1)
        plot!(tpath2[end:-1:2,2],CC2)
        savefig("../Results/BestWorst/CumConservBest$(i).png")
        # Make cumaltive entropy productions
        CΔS1 = zeros(length(ΔS1))
        CΔS2 = zeros(length(ΔS2))
        for i = 1:length(C1)
            CΔS1[i] = sum(ΔS1[1:i])
            CΔS2[i] = sum(ΔS2[end+1-i:end])
        end
        plot(tpath1[2:end,2],ΔS1)
        plot!(tpath2[end:-1:2,2],ΔS2[end:-1:1])
        savefig("../Results/BestWorst/EntPBest$(i).png")
        plot(tpath1[2:end,2],CΔS1)
        plot!(tpath2[end:-1:2,2],CΔS2)
        savefig("../Results/BestWorst/CumEntPBest$(i).png")
    end
    # Now read in worst trajectories and plot
    for i = 1:3#length(worst)
        infile1 = "../Results/Fig3Data/Traj/$(worst[i])$(ARGS[1])A2B.csv"
        infile2 = "../Results/Fig3Data/Traj/$(worst[i])$(ARGS[1])B2A.csv"
        # now should read in 1st path
        l = countlines(infile1)
        w = 2
        path1 = zeros(l,w)
        open(infile1, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 1
            for line in eachline(in_file)
                # parse line by finding commas
                L = length(line)
                comma = fill(0,w+1)
                m = 1
                for i = 1:L
                    if line[i] == ','
                        m += 1
                        comma[m] = i
                    end
                end
                comma[end] = L+1
                for i = 1:w
                    path1[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # now should read in 2nd path
        l = countlines(infile2)
        w = 2
        path2 = zeros(l,w)
        open(infile2, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            k = 1
            for line in eachline(in_file)
                # parse line by finding commas
                L = length(line)
                comma = fill(0,w+1)
                m = 1
                for i = 1:L
                    if line[i] == ','
                        m += 1
                        comma[m] = i
                    end
                end
                comma[end] = L+1
                for i = 1:w
                    path2[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                end
                k += 1
            end
        end
        # Calculate and then plot velocities
        v1, tpath1, f1, C1, ΔS1 = speed(path1,ps[worst[i],:])
        v2, tpath2, f2, C2, ΔS2 = speed(path2,ps[worst[i],:])
        plot(path1[:,2],path1[:,1])
        plot!(path2[:,2],path2[:,1])
        savefig("../Results/BestWorst/Worst$(i).png")
        # plot(tpath1[2:end,2],abs.(v1))
        # plot!(tpath2[2:end,2],abs.(v2))
        plot(tpath1[2:end,2],v1)
        plot!(tpath2[2:end,2],v2)
        savefig("../Results/BestWorst/VelocityWorst$(i).png")
        # plot(tpath1[2:end,2],abs.(f1))
        # plot!(tpath2[2:end,2],abs.(f2))
        plot(tpath1[2:end,2],f1)
        plot!(tpath2[2:end,2],f2)
        savefig("../Results/BestWorst/ForceWorst$(i).png")
        plot(tpath1[2:end,2],C1)
        plot!(tpath2[2:end,2],C2)
        savefig("../Results/BestWorst/ConservWorst$(i).png")
        # Make cumulative actions
        CC1 = zeros(length(C1))
        CC2 = zeros(length(C2))
        for i = 1:length(C1)
            CC1[i] = sum(C1[1:i])
            CC2[i] = sum(C2[end+1-i:end])
        end
        plot(tpath1[2:end,2],CC1)
        plot!(tpath2[end:-1:2,2],CC2)
        savefig("../Results/BestWorst/CumConservWorst$(i).png")
        # Make cumaltive entropy productions
        CΔS1 = zeros(length(ΔS1))
        CΔS2 = zeros(length(ΔS2))
        for i = 1:length(C1)
            CΔS1[i] = sum(ΔS1[1:i])
            CΔS2[i] = sum(ΔS2[end+1-i:end])
        end
        plot(tpath1[2:end,2],ΔS1)
        plot!(tpath2[end:-1:2,2],ΔS2[end:-1:1])
        savefig("../Results/BestWorst/EntPWorst$(i).png")
        plot(tpath1[2:end,2],CΔS1)
        plot!(tpath2[end:-1:2,2],CΔS2)
        savefig("../Results/BestWorst/CumEntPWorst$(i).png")
    end
    return(nothing)
end

@time main()
