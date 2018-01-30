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

# functions to calculate the nullclines for A at a given point x = (B, S, W)
function nullcline1(x)
    A1 = (1/(kmin + K))*(k*x[2]*(r/(r + f*x[1]^2)) + Kmin*x[3])
    return(A1)
end

function nullcline2(x)
    A2 = Complex((r/f)*((k*x[2]/((kmin + K)*x[1] - Kmin*x[3])) - 1))^0.5 # only valid if ϕ = q/k
    erA2 = 10.0^50 # Not sure if this is the best way to do this
    # using a large number in order to weight very heavily against it
    if imag(A2) != 0
        return(erA2)
        print("WARNING: In Dodgy Region (Nullcline2) Could Effect Validity Of Results\n")
    else
        return(real(A2))
    end
end

function nullcline3(x)
    A3 = (Kmin/K)*(1 + ϕ)*x[3] + (F/K) - ϕ*x[1] # only valid if ϕ = q/k
    return(A3)
end

function nullcline4(x)
    A4 = Complex((r/f)*((k*ϕ*x[2]/(F + ϕ*kmin*x[1] - k*x[2]*r/(r + f*x[1]^2)) - 1)))^0.5 # only valid if ϕ = q/k
    erA4 = 10.0^50
    if imag(A4) != 0
        return(erA4)
        print("WARNING: In Dodgy Region (Nullcline4) Could Effect Validity Of Results\n")
    else
        return(real(A4))
    end
end

# single function containing all the equations to minimize
function f!(F, x)
    F[1] = nullcline1(x[2:4]) - x[1]
    F[2] = nullcline2(x[2:4]) - x[1]
    F[3] = nullcline3(x[2:4]) - x[1]
    F[4] = nullcline4(x[2:4]) - x[1]
    F[5] = N - x[1] - x[2] - x[3] - x[4]
end

# jacobian should think a bit about how I do this
# x[1] = A, x[2] = B, x[3] = S, x[4] = W, x[5] = N
function j!(J, x)
    J[1, 1] = -1
    J[1, 2] = -2*f*x[2]*k*x[3]*r/((kmin + K)*(r + f*x[2]^2)^2)
    J[1, 3] = k*r/((kmin + K)*(r + f*x[2]^2))
    J[1, 4] = Kmin/(kmin + K)
    J[1, 5] = 0
    J[2, 1] = -1
    denom = (kmin + K)*x[2] - Kmin*x[4]
    # Need to think about how to deal with complex numbers here
    A = -((r/f)^0.5)*k*x[3]/(2*(denom)*(Complex(k*x[3]/(denom) - 1)^0.5))
    if imag(A) != 0
        print("Warning: 2nd Row Of Jacobian Found In Dodgy Region\n")
    end
    A = real(A)
    J[2, 2] = -(kmin + K)*A/(denom)
    J[2, 3] = A
    J[2, 4] = Kmin*A/(denom)
    J[2, 5] = 0
    J[3, 1] = -1
    J[3, 2] = -ϕ
    J[3, 3] = 0
    J[3, 4] = (Kmin/K)*(1 + ϕ)
    J[3, 5] = 0
    J[4, 1] = -1
    denom2 = F + ϕ*kmin*x[2] - r*k*x[3]/(r + f*x[2]^2)
    # And how to deal with complex numbers here
    B = ((r/f)^0.5)*ϕ*k/(2*(Complex(ϕ*k*x[3]/(denom2) - 1)^0.5))
    if imag(B) != 0
        print("Warning: 4th Row Of Jacobian Found In Dodgy Region\n")
    end
    B = real(B)
    J[4, 2] = B*x[3]*(ϕ*kmin + 2*r*f*x[2]*k*x[3]/((r + f*x[2]^2)^2))/(denom2^2)
    J[4, 3] = B*(1/(denom2) + k*r/((r + f*x[2]^2)*(denom2^2)))
    J[4, 4] = 0
    J[4, 5] = 0
    J[5, 1] = -1
    J[5, 2] = -1
    J[5, 3] = -1
    J[5, 4] = -1
    J[5, 5] = 0
end



function find_zeros()
    # system of 6 equations to find the minimum of
    xs = nlsolve(f!, j!, [200.0, 200.0, 800.0, 800.0, N], iterations = 1000)
    #xs = nlsolve(f!, [200.0, 200.0, 800.0, 800.0, N], iterations = 1000)
    print(xs)
    print("\n")
    n1 = nullcline1(xs.zero[2:4]) - xs.zero[1]
    n2 = nullcline2(xs.zero[2:4]) - xs.zero[1]
    n3 = nullcline3(xs.zero[2:4]) - xs.zero[1]
    n4 = nullcline4(xs.zero[2:4]) - xs.zero[1]
    print(n1+n2+n3+n4)
    print("\n")
end

@time find_zeros()
