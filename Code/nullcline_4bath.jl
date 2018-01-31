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
function nullclineS1(x)
    S1 = (r + f*x[2]^2)*((kmin + K)*x[1] - Kmin*x[3])/(k*r)
    return(S1)
end

function nullclineS2(x)
    S2 = (r + f*x[1]^2)*((kmin + K)*x[2] - Kmin*x[3])/(k*r)
    return(S2)
end

function nullclineS3(x)
    S3 = (kmin*(x[1] + ϕ*x[2]) + F)/(k*(r/(r + f*x[2]^2) + ϕ*r/(r + f*x[1]^2)))
    return(S3)
end

function nullclineW1(x)
    W1 = k*x[4]*r/(Kmin*(r + f*x[2]^2)) - (kmin + K)*x[1]/Kmin
    return(W1)
end

function nullclineW2(x)
    W2 = k*x[4]*r/(Kmin*(r + f*x[1]^2)) - (kmin + K)*x[2]/Kmin
    return(W2)
end

function nullclineW3(x)
    W3 = (K*(x[1] + ϕ*x[2]) - F)/(Kmin*(1 + ϕ))
    return(W3)
end

# single function containing all the equations to minimize
# x[1] = A, x[2] = B, x[3] = W, x[4] = S
# function termed F1 to avoid name space conflict with global constant F
function f!(F1, x)
    F1[1] = nullclineS1(x) - x[4]
    F1[2] = nullclineS2(x) - x[4]
    F1[3] = nullclineS3(x) - x[4]
    F1[4] = nullclineW1(x) - x[3]
    F1[5] = nullclineW2(x) - x[3]
    F1[6] = nullclineW3(x) - x[3]
    F1[7] = N - x[1] - x[2] - x[3] - x[4] # constraint term
end



function find_zeros()
    # system of 7 equations to find the minimum of
    xs = nlsolve(f!, [100.0, 700.0, 600.0, 600.0, 0.0, 0.0, 0.0], iterations = 1000) # three zeros are for unused variables
    # AiBi = N/10
    # Bi = AiBi/10
    # Ai = AiBi - Bi
    # Si = (N - AiBi)/10
    # Wi = 9*(N - AiBi)/10
    # xs = nlsolve(f!, [Ai, Bi, Si, Wi, N])
    print(xs)
    print("\n")
    n1 = nullclineS1(xs.zero) - xs.zero[4]
    n2 = nullclineS2(xs.zero) - xs.zero[4]
    n3 = nullclineS3(xs.zero) - xs.zero[4]
    n4 = nullclineW1(xs.zero) - xs.zero[3]
    n5 = nullclineW2(xs.zero) - xs.zero[3]
    n6 = nullclineW3(xs.zero) - xs.zero[3]
    n7 = N - xs.zero[1] - xs.zero[2] - xs.zero[3] - xs.zero[4]
    print(n1)
    print("\n")
    print(n2)
    print("\n")
    print(n3)
    print("\n")
    print(n4)
    print("\n")
    print(n5)
    print("\n")
    print(n6)
    print("\n")
    print(n7)
    print("\n")
    print(n1 + n2 + n3 + n4 + n5 + n6 + n7)
    print("\n")
end

@time find_zeros()
