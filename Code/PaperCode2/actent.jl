#!/usr/bin/env julia
# actent.jl
# A script to read in the parameters and paths and then use to determine
# entropy produced along the path and action of forward and reverse paths
#
# Author: Jacob Cook
# Date: December 2018

using SymEngine

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

# Function to calculate the action of a given path
function Ŝ(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},NG::Int64)
    S = zeros(NG-1)
    S[1] = (3/(2*NG))*(xprim[2,1]*ϑs[2,1] + xprim[2,2]*ϑs[2,2])
    # Not excluding the midpoint as the contribution is vanishing
    # Might have to rethink this for the 4 species case
    for i = 3:NG-1
        S[i-1] += (1/NG)*(xprim[i,1]*ϑs[i,1] + xprim[i,2]*ϑs[i,2])
    end
    S[NG-1] += (3/(2*NG))*(xprim[NG,1]*ϑs[NG,1] + xprim[NG,2]*ϑs[NG,2])
    return(S)
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

# function to take time discretised path and return Action and total entropy production based on it
function act(x::Array{Float64,2},Tp::Float64,b::Array{SymEngine.Basic,1},Dmin::Array{SymEngine.Basic,2})
    # find time step etc
    N = size(x,1)-1
    dt = Tp/(N)
    Ac = 0
    ΔS = 0
    A, B, y1, y2 = symbols("A B y1 y2")
    # temporary bs and Dmins to store values
    bt = Array{Float64,1}(undef,2)
    Dm = Array{Float64,2}(undef,2,2)
    for i = 1:N
        # need to find velocity for each point
        qdot = (x[i+1,:] .- x[i,:])/(dt)
        # and need to find midpoint
        posA = (x[i+1,1] + x[i,1])/(2)
        posB = (x[i+1,2] + x[i,2])/(2)
        # now substitue into b and Dmin
        for j = 1:2
            bt[j] = subs(b[j],A=>posA,B=>posB) |> float
            for k = 1:2
                Dm[j,k] = subs(Dmin[j,k],A=>posA,B=>posB,y1=>qdot[1],y2=>qdot[2]) |> float
            end
        end
        # now calculate the update to A and ΔS from this point
        for j = 1:2
            for k = 1:2
                Ac += 0.5*(qdot[j] - bt[j])*Dm[j,k]*(qdot[k] - bt[k])*dt
                ΔS += 2*qdot[j]*Dm[j,k]*bt[k]*dt
            end
        end
    end
    return(Ac,ΔS)
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
    # 2 for loops to loop over every file
    for i = 1:l
        for j = 1:2
            if j == 1
                infile = "../Results/Fig3Data/Traj2/$(i)$(ARGS[1])A2B.csv"
                outfile = "../Results/Fig3Data/Traj2/$(i)$(ARGS[1])A2BD.csv"
            else
                infile = "../Results/Fig3Data/Traj2/$(i)$(ARGS[1])B2A.csv"
                outfile = "../Results/Fig3Data/Traj2/$(i)$(ARGS[1])B2AD.csv"
            end
            # check if file
            if isfile(infile)
                # now should read in path
                l = countlines(infile)
                w = 2
                path = zeros(l,w)
                open(infile, "r") do in_file
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
                            path[k,i] = parse(Float64,line[(comma[i]+1):(comma[i+1]-1)])
                        end
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
                # now use a function that takes the time discretised path and
                # finds the action in a more conventional manner and then can also get entropy production from this
                Act, ΔS = act(path2,Tp,b,Dmin)
                # and write out the new data to a new file
                out_file = open(outfile, "w")
                line = "$(Tp),$(ActS),$(Act),$(ΔS)\n"
                write(out_file, line)
                close(out_file)
            else # Tell users about missing files
                if j == 1
                    println("No A to B path for $(i)")
                else
                    println("No B to A path for $(i)")
                end
            end
        end
    end
    return(nothing)
end

@time main()
