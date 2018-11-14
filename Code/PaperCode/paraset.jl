#!/usr/bin/env julia
# paraset.jl
# A script to generate and then output a set of parameters for 2 species gene switch
# these parameters are then output as a .csv file
#
# Author: Jacob Cook
# Date: November 2018

using Roots

# function to generate unifrom random numbers between from and to
function myrand(to::Float64,from::Float64)
    r = rand()*(to-from) + from
    return(r)
end

# Need to think of a good parameter randomisation scheme
function paras()
    works = false
    paras = zeros(10)
    while works == false
        # k randomly selected between 1.0 and 100.0
        k = myrand(100.0,1.0)
        kmin = k*myrand(0.1,0.001)
        # q then smaller than k by a ratio between 1.0 and 0.01
        q = k*myrand(1.0,0.01)
        qmin = q*myrand(0.1,0.001)
        # at this point I need to think more carefully to make sure these are parameters with steady states
        # K randomly selected between 0.1 and 10.0
        K = myrand(10.0,0.1)
        Kmin = K*myrand(0.1,0.001)
        # Q then smaller than K by a ratio between 1.0 and 0.01
        Q = K*myrand(1.0,0.01)
        Qmin = Q*myrand(0.1,0.001)
        r = 10000.0*rand()
        f = 10000.0*rand()
        # stopping condition
        # now gather parameters
        paras = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
        steady = steady3(paras)
        if steady == true
            works = true
        end
    end
    return(paras)
end

# function to test if parameter set allows three +ve stationary points
function steady3(ps)
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    # A1(x) = sqrt((r/f)*(q/((qmin+Q)*B - Qmin) - 1))
    A1(x) = real(sqrt(complex((ps[9]/ps[10])*(ps[3]/((ps[4]+ps[7])*x - ps[8]) - 1))))
    # A2(x) = (1/(kmin+K))*(kr/(r+f*x^2) + Kmin)
    A2(x) = (1/(ps[2]+ps[5]))*((ps[1]*ps[9])/(ps[9]+ps[10]*x^2) + ps[6])
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    q = ps[3]
    Q = ps[7]
    for i = 1:2
        bs1 = fzeros(g, 0.0, 0.1) # catches zeros about the origin
        bs2 = fzeros(g, 0.1, 2*q/Q) # the upper bound here is potentially problematic
        bs = vcat(bs1,bs2)
        # bs = fzeros(g, 0, 15.0)
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
            return(true)
        end
    end
    return(false)
end

function main()
    println("Compiled, Starting script.")
    # First check that an argument for naming has been provided
    if length(ARGS) == 0
        println("Error: Need to provide an argument to name output with.")
        return(nothing)
    # then check that a system volume has been provided
    elseif length(ARGS) == 1
        println("Error: Need to provide an argument to set number of parameter sets.")
        return(nothing)
    end
    # Then take system volume Ω, check if provided value is integer
    N = 0
    try N = parse(Int64,ARGS[2])
    catch y
        if isa(y, ArgumentError) # would only really expect an argument error
            println("Error: Number of parameter sets has to be integer.")
            return(nothing)
        end
    end
    # ps = [ k, kmin, q, qmin, K, Kmin, Q, Qmin, r, f]
    ps = zeros(N,10) # matrix to store full parameter set
    for i = 1:N
        ps[i,:] = paras()
    end
    # Finally write out to file
    output_file = "../Results/Fig3Data/$(ARGS[1])para.csv"
    out_file = open(output_file, "w")
    for i = 1:size(ps,1)
        line = ""
        for j = 1:size(ps,2)
            line *= "$(ps[i,j]),"
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
