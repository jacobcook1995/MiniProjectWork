#!/usr/bin/env julia
# Entp.jl
# script to read in a long trajectory and then calculate the entropy production of it

using Plots
using SymEngine
import GR

function readin()
    filename = "../Results/$(ARGS[1]).csv"
    # count lines
    no_lins = countlines(filename)
    # then preallocate arrays
    vars = zeros(4,no_lins)
    times = zeros(no_lins)
    open(filename, "r") do in_file
        k = 0
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            k += 1
            # parse line by finding commas
            comma = fill(0,4)
            j = 0
            L = length(line)
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            vars[1,k] = parse(Float64, line[1:(comma[1] - 1)])
            vars[2,k] = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
            vars[3,k] = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
            vars[4,k] = parse(Float64, line[(comma[3] + 1):(comma[4] - 1)])
            times[k] = parse(Float64, line[(comma[4] + 1):L])
        end
    end
    return(vars,times)
end
# Vector of functions from MAP case
function f!(G, x, ps)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin = symbols("K k Q q kmin qmin f r F Kmin Qmin")
    sym = Array{SymEngine.Basic}(4)
    sym[1] = k*x[3]*r/(r+f*x[2]^2) - K*x[1] - kmin*x[1] + Kmin*x[4]
    sym[2] = q*x[3]*r/(r+f*x[1]^2) - Q*x[2] - qmin*x[2] + Qmin*x[4]
    sym[3] = -k*x[3]*r/(r + f*x[2]^2) - q*x[3]*r/(r + f*x[1]^2) + kmin*x[1] + qmin*x[2] + F
    sym[4] = K*x[1] + Q*x[2] - Kmin*x[4] - Qmin*x[4] - F
    for i = 1:4
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        G[i] = subs(sym[i], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11]) |> float
    end
    return G
end

# Diffusion matrix from MAP case
function D!(D, x, ps)
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(4,4)
    e[1,2:4] = e[2,1] = e[2,3:4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    e[4,1] = -sqrt(K*A + Kmin*W)
    e[4,2] = -sqrt(Q*B + Qmin*W)
    e[4,4] = sqrt(F) # this is not as previously discussed, could be the error
    # Now do the transformations required
    eT = transpose(e)
    D = e*eT
    D = inv(D)
    for j = 1:4
        for i = 1:4
           D[i,j] = subs(D[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
           D[i,j] = subs(D[i,j], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11], A=>x[1], B=>x[2], S=>x[3], W=>x[4]) |> float
       end
    end
    return D
end

# Function to calculate the action for a given point on the path
function EntProd(vars::Array{Float64,2},τ::Array{Float64,1},Ω)
    # parameters hard coded into this function
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0 # Promoter switching
    r = 10.0
    F = 10.0
    Kmin = 10.0^-20
    Qmin = 10.0^-20
    ps = [ K; k; Q; q; kmin; qmin; f; r; F; Kmin; Qmin]
    termsthi = zeros(4,4)
    termsthif = zeros(4,4)
    termsf = zeros(4,4)
    Fqs = zeros(4)
    Ffs = zeros(4)
    Fs = 0
    Acts = 0
    Ents = 0
    # maybe can avoid the reallocation trap here
    ############################################################################
    h = [ 0.0; 0.0; 0.0; 0.0 ]
    thiv = [ 0.0; 0.0; 0.0; 0.0 ]
    d = [ 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 ]
    ############################################################################
    deltat = τ[2] - τ[1]
    posA = (vars[1,1] + vars[1,2])/2
    posB = (vars[2,1] + vars[2,2])/2
    posS = (vars[3,1] + vars[3,2])/2
    posW = (vars[4,1] + vars[4,2])/2
    h = f!(h, [posA posB posS posW], ps)
    h[3] -= F
    h[4] += F
    d = D!(d, [posA posB posS posW], ps)
    for j = 1:4
        thiv[j] = (vars[j,2] - vars[j,1])/(deltat)
    end
    # Now calculation step
    for k = 1:4
        for j = 1:4
            termsthi[j,k] = thiv[j]*d[j,k]*thiv[k]
            termsthif[j,k] = h[j]*d[j,k]*thiv[k]
            termsf[j,k] = h[j]*d[j,k]*h[k]
        end
        Fqs[k] = thiv[k]*d[k,3]*F - thiv[k]*d[k,4]*F
        Ffs[k] = h[k]*d[k,3]*F - h[k]*d[k,4]*F
    end
    Fs = F*d[3,3]*F + F*d[4,4]*F

    # Now use these terms to calculate overall action along path and the entropy production
    Acts = Fs*(deltat/2)
    for k = 1:4
        for j = 1:4
            Acts += termsthi[j,k]*(deltat/2)
            Acts -= termsthif[j,k]*deltat
            Acts += termsf[j,k]*(deltat/2)
            Ents += termsthif[j,k]*(2*deltat)
        end
        Acts -= Fqs[k]*deltat
        Acts += Ffs[k]*deltat
        Ents += Fqs[k]*(2*deltat) # comment this out to see the effect at somepoint
    end
    return(Acts,Ents,deltat,thiv,Ffs,Fqs,Fs,termsthi,termsf,termsthif)
end

function main()
    # check if argument is provided
    if length(ARGS) >= 1
        vars, times = readin()
    else
        println("Need to provide an argument")
    end
    # Now need to rescale the trajectories
    Ω = 60
    vars = vars/Ω
    # Now need to do the calculation of entropy production
    At = St = 0
    A = zeros(4999999)
    S = zeros(4999999)
    for i = 2:5000000#length(times)
        A[i-1], S[i-1] = EntProd(vars[:,i-1:i],times[i-1:i],Ω)
        At += A[i-1]
        St += S[i-1]
        if i % 50000 == 0
            println("$(i),$(At/times[i-1]),$(St/times[i-1])")
        end
    end
    plot(A)
    savefig("../Results/A$(ARGS[1]).png")
    plot(S)
    savefig("../Results/S$(ARGS[1]).png")
    return(nothing)
end

@time main()
