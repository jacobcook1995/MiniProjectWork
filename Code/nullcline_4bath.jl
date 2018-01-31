#!/usr/bin/env julia
# nullcline_4bath.jl
# A script to find the nullclines and thus the steady states in the 4 bath case
# Hopefully this script will visulize them and the flow, though this will be hard
# as 4 dimensions are being considered

# Add package for plotting
using Plots
using NLsolve
import GR

# Parameters
const Ω = 30 # system size, less important parameter now but should still retain
const ϕ = 0.1 # ratio ϕ = q/k
const k = 100 # steady state for A=k/K=1
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const K = k/Ω # K=k'
const Kmin = 10.0^-20
const q = k*ϕ # steady state for B=q/Q=1
const qmin = 10.0^-20
const Q = q/Ω # Q=q'
const Qmin = 10.0^-20
const f = 10 # Promoter switching
const r = 10
const F = 250 # removal rate
const N = 2000.0 # density within the system

# functions to calculate the nullclines for A at a given point
# x[1] = A, x[2] = B, x[3] = W, x[4] = S
function nullclineA(x)
    A = (k*x[4]*r/(r + f*x[2]^2) + Kmin*x[3])/(kmin + K)
    return(A)
end

function nullclineB(x)
    B = (k*x[4]*r/(r + f*x[1]^2) + Kmin*x[3])/(kmin + K)
    return(B)
end

function Wdot(x)
    wdot = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3] - F
    return(wdot)
end

function nullclineS(x)
    S = (kmin*(x[1] + ϕ*x[2]) + F)/(k*(r/(r + f*x[2]^2) + ϕ*r/(r + f*x[1]^2)))
    return(S)
end

# single function containing all the equations to minimize
# x[1] = A, x[2] = B, x[3] = W, x[4] = S
# function termed F1 to avoid name space conflict with global constant F
function f!(F1, x)
    F1[1] = nullclineA(x) - x[1]
    F1[2] = nullclineB(x) - x[2]
    F1[3] = nullclineS(x) - x[4]
    F1[4] = x[1] + x[2] + x[3] + x[4] - N # constraint term
end


# function takes in initial values and either throws an error or returns the solution it finds
function find_zeros(Ai, Bi, Si, Wi)
    xs = nlsolve(f!, [Ai, Bi, Si, Wi])
    nW = Wdot(xs.zero)
    if abs(nW) <= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == false
        return(xs, nW)
    elseif abs(nW) >= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == false
        throw(DomainError) # just about reasonable
        error("PointChoiceError: not a valid initial point, w not valid")
    elseif abs(nW) <= 10.0^-10 && converged(xs) == false && any(x -> x < 0.0, xs.zero) == false
        throw(BoundsError)
        error("PointChoiceError: not a valid initial point cannot converge")
    elseif abs(nW) <= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == true
        throw(InexactError) # dubious use of exceptions but I can't be bothered to make my own at the minute
        error("PointChoiceError: not a valid initial point result has negative components")
    else
        throw(TypeError)
        error("PointChoiceError: not a valid initial point fails multiple conditions")
    end
end

# function to search initial posistion space to find any other points that work bar the steady state
function point_search()
    # High A
    # Ai = N/10
    # Bi = N/10
    # Si = 2*N/5
    # Wi = 2*N/5
    # Higher B
    #Ai = N/200
    #Bi = 2*N/5
    #Si = 119*N/400
    #Wi = 119*N/400
    for i = 1:10
        rA = rand()
        rB = rand()
        rS = rand()
        rW = rand()
        rs = [ rA; rB; rS; rW ]
        rs = rs/sum(rs)
        Ai = N*rs[1]
        Bi = N*rs[2]
        Si = N*rs[3]
        Wi = N*rs[4]
        try find_zeros(Ai, Bi, Si, Wi)
        catch y
            if isa(y, DomainError) # the errors here don't really mean anything I'm just too lazy to make my own
                print("Does not work due to waste not being steady\n")
            elseif isa(y, BoundsError)
                print("Does not work due to not convergance\n")
            elseif isa(y, InexactError)
                print("Does not work due to negativity\n")
            elseif isa(y, TypeError)
                print("Does not work due to everything\n")
            end
        finally
            print("a\n")
        end
    end
end

@time point_search()
