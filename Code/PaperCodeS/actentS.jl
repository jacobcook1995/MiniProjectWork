#!/usr/bin/env julia
# actentS.jl
# A script to read in the parameters and paths and then use to determine
# entropy produced along the path and action of forward and reverse paths
# This is now done for the 1 species Schlögl model
#
# Author: Jacob Cook
# Date: December 2018

using SymEngine
using LaTeXStrings
using Plots
import PyPlot

# make a symbolic diffusion matrix
function Ds()
    X, k1, K1, k2, K2 = symbols("X k1 K1 k2 K2")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,1,1)
    e[1,1] = sqrt(k2*X^3 + K2*X^2 + K1*X + k1)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    Dmin[1,1] = expand(Dmin[1,1])
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    X, k1, K1, k2, K2 = symbols("X k1 K1 k2 K2")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,1)
    b[1] = -k2*X^3 + K2*X^2 - K1*X + k1
    return(b)
end

# vector of in forces
function f1!(f::Float64,x::Float64,ps::Array{Float64,1})
    f = ps[1] - ps[2]*x
    return(f)
end

# vector of out forces
function f2!(f::Float64,x::Float64,ps::Array{Float64,1})
    f = ps[3]*x^3 - ps[4]*x^2
    return(f)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y = symbols("y")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1]
    den = 0
    den += y*Dmin[1,1]*y
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

# function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y = symbols("y")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,1)
    c[1] = λ*y - b[1]
    ϑ = Array{SymEngine.Basic,1}(undef,1)
    ϑ[1] = Dmin[1,1]*c[1]
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
    k1, K1, k2, K2 = symbols("k1 K1 k2 K2")
    # now perform substitutions
    # ps = [ k1, K1, k2, K2 ]
    λ = subs(λ, k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    ϑ[1] = subs(ϑ[1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    b[1] = subs(b[1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    Dmin[1,1] = subs(Dmin[1,1], k1=>ps[1], K1=>ps[2], k2=>ps[3], K2=>ps[4])
    return(ϑ,λ,b,Dmin)
end

# function to generate the variables needed for a given path
function genvars(x::Array{Float64,1},λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int64,Nmid::Int64)
    # define neccesary symbols
    X, y = symbols("X y")
    # calculate velocities
    xprim = fill(NaN, NG+1)
    for i = 2:NG
        xprim[i] = (x[i+1] - x[i-1])/(2/NG)
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λs[i] = subs(λt, X=>x[i], y=>xprim[i]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λs[i] = subs(λt, X=>x[i], y=>xprim[i]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1)
    ϑt = Array{SymEngine.Basic,1}(undef,1)
    for i = 2:NG
        ϑt[1] = subs(ϑ[1], X=>x[i])
        ϑs[i] = subs(ϑt[1], y=>xprim[i]) |> float
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# Function to calculate the action of a given path
function Ŝ(x::Array{Float64,1},xprim::Array{Float64,1},λs::Array{Float64,1},ϑs::Array{Float64,1},λprim::Array{Float64,1},NG::Int64)
    S = zeros(NG-1)
    S[1] = (3/(2*NG))*(xprim[2]*ϑs[2])
    # Not excluding the midpoint as the contribution is vanishing
    # Might have to rethink this for the 4 species case
    for i = 3:NG-1
        S[i-1] += (1/NG)*(xprim[i]*ϑs[i])
    end
    S[NG-1] += (3/(2*NG))*(xprim[NG]*ϑs[NG])
    return(S)
end

# function to find the times of each point
function times(x::Array{Float64,1},xprim::Array{Float64,1},λs::Array{Float64,1},ϑs::Array{Float64,1},λprim::Array{Float64,1},NG::Int64)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::Array{Float64,1},x::Array{Float64,1},NG::Int64,NM::Int64)
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
    path = zeros(NM+1)
    path[1] = x[1]
    path[NM+1] = x[NG+1]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        x₀ = x[one]
        x₁ = x[two]
        path[i] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
    end
    return(path)
end

# function to take time discretised path and return Action and total entropy production based on it
function act(x::Array{Float64,1},Tp::Float64,b::Array{SymEngine.Basic,1},Dmin::Array{SymEngine.Basic,2},ps::Array{Float64,1},Ns::Int64)
    # find time step etc
    N = size(x,1)-1
    dt = Tp/(N)
    Acf = zeros(N)
    ΔSf = zeros(N)
    KE = zeros(N)
    PE = zeros(N)
    prods = zeros(N)
    flows = zeros(N)
    fs = zeros(N)
    qs = zeros(N)
    X, y = symbols("X y")
    desd = false
    f1 = 0.0
    f2 = 0.0
    # temporary bs and Dmins to store values
    bt = Array{Float64,1}(undef,1)
    Dm = Array{Float64,2}(undef,1,1)
    for i = 1:N
        # condition to switch
        if i == Ns
            desd = true
        end
        # need to find velocity for each point
        qdot = (x[i+1] - x[i])/(dt)
        # and need to find midpoint
        posX = (x[i+1] + x[i])/(2)
        # now substitue into b and Dmin
        bt[1] = subs(b[1],X=>posX) |> float
        Dm[1,1] = subs(Dmin[1,1],X=>posX,y=>qdot) |> float
        f1 = f1!(f1,posX,ps)
        f2 = f2!(f2,posX,ps)
        # now calculate the update to A and ΔS from this point
        Acf[i] += 0.5*(qdot - bt[1])*Dm[1,1]*(qdot - bt[1])*dt
        ΔSf[i] += 2*qdot*Dm[1,1]*bt[1]*dt
        KE[i] += 0.5*qdot*Dm[1,1]*qdot*dt
        PE[i] += 0.5*bt[1]*Dm[1,1]*bt[1]*dt
        fs[i] += bt[1]
        qs[i] += qdot
        if desd == false
            prods[i] += 4*f1*Dm[1,1]*f2*dt
            flows[i] += 2*(f1*Dm[1,1]*f1 + f2*Dm[1,1]*f2)*dt
        else
            prods[i] += 2*(f1*Dm[1,1]*f1 + f2*Dm[1,1]*f2)*dt
            flows[i] += 4*f1*Dm[1,1]*f2*dt
        end
    end
    Ac = sum(Acf)
    ΔS = sum(ΔSf)
    return(Ac,ΔS,Acf,ΔSf,KE,PE,prods,flows,fs,qs)
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
    # 2 for loops to loop over every file
    for i = 1:l
        for j = 1:2
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2l.csv"
                outfile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2lD.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2h.csv"
                outfile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2hD.csv"
            end
            # check if file
            if isfile(infile)
                # now should read in path
                l = countlines(infile)
                w = 1
                path = zeros(l)
                open(infile, "r") do in_file
                    # Use a for loop to process the rows in the input file one-by-one
                    k = 1
                    for line in eachline(in_file)
                        # just parse entire line
                        path[k] = parse(Float64,line[1:end])
                        k += 1
                    end
                end
                # generate symbolic objects
                ϑ, λ, b, Dmin = gensyms(ps[i,:])
                # define NG and Nmid and use to find variables
                NG = NM = l - 1
                Nmid = convert(Int64, ceil((NG+1)/2))
                x, xprim, λs, ϑs, λprim = genvars(path,λ,ϑ,NG,Nmid)
                # use function Ŝ to find the action associated with this path
                λs[1] = λs[2]
                λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
                λs[end] = λs[end-1]
                # find and save action from geometric method
                S = Ŝ(x,xprim,λs,ϑs,λprim,NG)
                ActS = sum(S)
                tims2 = times(x,xprim,λs,ϑs,λprim,NG)
                # save time of path
                Tp = tims2[end]
                path2 = timdis(tims2,x,NG,NM)
                Ns = convert(Int64,NM/2)
                # now use a function that takes the time discretised path and
                # finds the action in a more conventional manner and then can also get entropy production from this
                Act, ΔS = act(path2,Tp,b,Dmin,ps[i,:],Ns)
                # and write out the new data to a new file
                out_file = open(outfile, "w")
                line = "$(Tp),$(ActS),$(Act),$(ΔS)\n"
                write(out_file, line)
                close(out_file)
            else # Tell users about missing files
                if j == 1
                    println("No high to low path for $(i)")
                else
                    println("No low to high path for $(i)")
                end
            end
        end
    end
    return(nothing)
end

function plotting()
    println("Compiled, Starting plotting script.")
    flush(stdout)
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    elseif length(ARGS) == 1
        println("Error: Need to provide an Integer to choose a file with.")
        return(nothing)
    end
    # Check integer has been provided
    N = 0
    int = true
    for i = 1:length(ARGS[2])
        if ~isnumeric(ARGS[2][i])
            int = false
            break
        end
    end
    if int == true
         N = parse(Int64,ARGS[2])
    else
        println("Error: Must provide an integer file number.")
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
    if N > l
        println("Error: File number provided is too large")
        return(nothing)
    end
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
    # 2 for loops to loop over every file
    pyplot()
    p1 = plot(dpi=300,title="Schlögl Paths",xlabel="Time (t)",ylabel="Concentration x")
    plot!(p1,titlefontsize=20,guidefontsize=16,legendfontsize=12,tickfontsize=14)
    mag1 = 10.0^3
    Smag1 = L"(10^{-3}\,1/\Omega)"
    p2 = plot(dpi=300,title="Schlögl Action Contributions",xlabel="Concentration x")
    plot!(p2,titlefontsize=20,guidefontsize=16,legendfontsize=12,ylabel="Action Contributions $(Smag1)",tickfontsize=14)
    p3 = plot(dpi=300,title="Schlögl Production and Flow Terms",xlabel="Concentration x")
    qlab = L"\dot{q}"
    plot!(p3,titlefontsize=20,guidefontsize=16,legendfontsize=12,ylabel="Action Contributions $(Smag1)",tickfontsize=14)
    p4 = plot(dpi=300,title="Comparison of f and $(qlab) for Schlögl model",xlabel="Concentration x",ylabel="Magnitude (Copy Number/s)")
    p5 = plot(dpi=300,grid=false,ticks=false,xaxis=false,yaxis=false)
    for i = N
        for j = 1:2
            if j == 1
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])h2l.csv"
            else
                infile = "../Results/Fig3DataS/Traj/$(i)$(ARGS[1])l2h.csv"
            end
            # check if file
            if isfile(infile)
                # now should read in path
                l = countlines(infile)
                w = 1
                path = zeros(l)
                open(infile, "r") do in_file
                    # Use a for loop to process the rows in the input file one-by-one
                    k = 1
                    for line in eachline(in_file)
                        # just parse entire line
                        path[k] = parse(Float64,line[1:end])
                        k += 1
                    end
                end
                # generate symbolic objects
                ϑ, λ, b, Dmin = gensyms(ps[i,:])
                # define NG and Nmid and use to find variables
                NG = NM = l - 1
                Nmid = convert(Int64, ceil((NG+1)/2))
                x, xprim, λs, ϑs, λprim = genvars(path,λ,ϑ,NG,Nmid)
                # use function Ŝ to find the action associated with this path
                λs[1] = λs[2]
                λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
                λs[end] = λs[end-1]
                # find and save action from geometric method
                S = Ŝ(x,xprim,λs,ϑs,λprim,NG)
                ActS = sum(S)
                tims2 = times(x,xprim,λs,ϑs,λprim,NG)
                # save time of path
                Tp = tims2[end]
                path2 = timdis(tims2,x,NG,NM)
                if j == 1
                    plot!(p1,tims2,path2,label="High → Low",color=1)
                    # now want to add arrows along the path indicating directions
                    plot!(p1,tims2[100:101],path2[100:101],arrow=0.4,color=1,label="")
                    plot!(p1,tims2[450:451],path2[450:451],arrow=0.4,color=1,label="")
                else
                    plot!(p1,tims2,path2,label="Low → High",color=2)
                    plot!(p1,tims2[150:151],path2[150:151],arrow=0.4,color=2,label="")
                    plot!(p1,tims2[500:501],path2[500:501],arrow=0.4,color=2,label="")
                end
                # find point that path reaches saddle Ns
                Ns = 0
                if path2[1] == steads[N,1]
                    for i = 1:length(path2)
                        if path2[i] <= steads[N,2]
                            Ns = i
                            break
                        end
                    end
                elseif path2[1] == steads[N,3]
                    for i = 1:length(path2)
                        if path2[i] >= steads[N,2]
                            Ns = i
                            break
                        end
                    end
                else
                    println("Error: Should start at sone of the two points")
                    error()
                end
                # now use a function that takes the time discretised path and
                # finds the action in a more conventional manner and then can also get entropy production from this
                Act, ΔS, Af, ΔSf, KE, PE, prods, flows, fs, qs = act(path2,Tp,b,Dmin,ps[i,:],Ns)
                if j == 1
                    plot!(p2,path2[1:end-1],Af*mag1,label="",color=1,style=:dash)
                    plot!(p2,path2[1:end-1],ΔSf*mag1,label="",color=2,style=:dash)
                    plot!(p2,path2[1:end-1],KE*mag1,label="",color=3,style=:dash)
                    plot!(p2,path2[1:end-1],PE*mag1,label="",color=4,style=:dash)
                    plot!(p3,path2[1:end-1],prods*mag1,label="",color=1,style=:dash)
                    plot!(p3,path2[1:end-1],flows*mag1,label="",color=2,style=:dash)
                    plot!(p3,path2[1:end-1],(prods.-flows)*mag1,label="",color=3,style=:dash)
                    plot!(p4,path2[1:end-1],fs,label=L"f",color=1)
                    plot!(p4,path2[1:end-1],qs,label=L"\dot{q}",color=2)
                    plot!(p5,path2[1:end-1],fs,label="",color=1)
                    plot!(p5,path2[1:end-1],qs,label="",color=2)
                    annotate!(p5,[(path2[68],fs[68]-0.5,text(L"f",20)),(path2[68],qs[68]+0.5,text(L"\dot{q}",20))])
                else
                    plot!(p2,path2[1:end-1],Af*mag1,label=L"\mathcal{A}",color=1)
                    plot!(p2,path2[1:end-1],ΔSf*mag1,label=L"\Delta S",color=2)
                    plot!(p2,path2[1:end-1],KE*mag1,label="KE",color=3)
                    plot!(p2,path2[1:end-1],PE*mag1,label="PE",color=4)
                    plot!(p3,path2[1:end-1],prods*mag1,label="EP",color=1)
                    plot!(p3,path2[1:end-1],flows*mag1,label="EF",color=2)
                    plot!(p3,path2[1:end-1],(prods.-flows)*mag1,label="EP-EF",color=3)
                end
            else # Tell users about missing files
                if j == 1
                    println("No high to low path for $(i)")
                else
                    println("No low to high path for $(i)")
                end
            end
        end
    end
    # add lines and then save
    hline!(p1,[steads[N,3]],color=:black,linestyle=:dash,label="")
    hline!(p1,[steads[N,2]],color=:black,linestyle=:dot,label="Saddle")
    hline!(p1,[steads[N,1]],color=:black,linestyle=:solid,label="")
    annotate!(p1,[(2.5,steads[N,3]+0.2,text("A",20)),(2.5,steads[N,1]-0.2,text("B",20))])
    savefig(p1,"../Results/Fig2Graphs/SchPath$(ARGS[2]).png")
    scatter!(p2,[steads[N,1]],[0.0],markersize=6,markercolor=:black,label="")
    scatter!(p2,[steads[N,2]],[0.0],markersize=5,markercolor=:black,markershape=:x,label="")
    scatter!(p2,[steads[N,3]],[0.0],markersize=6,markercolor=:white,label="")
    plot!(p2,ylims=(-8.0,8.0))
    savefig(p2,"../Results/Fig2Graphs/SchAct$(ARGS[2]).png")
    scatter!(p3,[steads[N,1]],[0.0],markersize=6,markercolor=:black,label="")
    scatter!(p3,[steads[N,2]],[0.0],markersize=5,markercolor=:black,markershape=:x,label="")
    scatter!(p3,[steads[N,3]],[0.0],markersize=6,markercolor=:white,label="")
    savefig(p3,"../Results/Fig2Graphs/SchfTerms$(ARGS[2]).png")
    scatter!(p4,[steads[N,1]],[0.0],markersize=6,markercolor=:black,label="")
    scatter!(p4,[steads[N,2]],[0.0],markersize=5,markercolor=:black,markershape=:x,label="")
    scatter!(p4,[steads[N,3]],[0.0],markersize=6,markercolor=:white,label="")
    savefig(p4,"../Results/Fig2Graphs/Schfqdot$(ARGS[2]).png")
    scatter!(p5,[steads[N,1]],[0.0],markersize=6,markercolor=:black,label="")
    scatter!(p5,[steads[N,2]],[0.0],markersize=5,markercolor=:black,markershape=:x,label="")
    scatter!(p5,[steads[N,3]],[0.0],markersize=6,markercolor=:white,label="")
    savefig(p5,"../Results/Fig2Graphs/Schredfqdot$(ARGS[2]).png")
    return(nothing)
end

# @time main()
@time plotting()
