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
    print(xs)
    print("\n")
    nW = Wdot(xs.zero)
    if abs(nW) <= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == false
        return(xs.zero)
    elseif abs(nW) >= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == false
        error("PointChoiceError: not a valid initial point, w not valid")
    elseif abs(nW) <= 10.0^-10 && converged(xs) == false && any(x -> x < 0.0, xs.zero) == false
        error("PointChoiceError: not a valid initial point cannot converge")
    elseif abs(nW) <= 10.0^-10 && converged(xs) == true && any(x -> x < 0.0, xs.zero) == true
        error("PointChoiceError: not a valid initial point result has negative components")
    else
        error("PointChoiceError: not a valid initial point fails all conditions")
    end
end

# High A
# Ai = N/10
# Bi = N/10
# Si = 2*N/5
# Wi = 2*N/5
# Higher B
Ai = N/200
Bi = 2*N/5
Si = 119*N/400
Wi = 119*N/400

@time X = find_zeros(Ai, Bi, Si, Wi)
print(X)
print("\n")
