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
const N = 2000 # number of elements in the system


# functions to calculate the nullclines for A at a given point x = (B, S, W)
function nullcline1(x)
    A1 = (1/(kmin + K))*(k*x[2]*(r/(r + f*x[1]^2)) + Kmin*x[3])
    return(A1)
end

function nullcline2(x)
    A2 = Complex((r/f)*((k*x[2]/((kmin + K)*x[1] - Kmin*x[3])) - 1))^0.5 # only valid if ϕ = q/k
    erA2 = Inf # Not sure if this is the best way to do this
    if imag(A2) != 0
        return(erA2)
    else
        return(real(A2))
    end
end

function nullcline3(x)
    A3 = (Kmin/K)*(1 + ϕ)*x[3] + (F/K) - ϕ*x[1] # only valid if ϕ = q/k
    return(A3)
end

function nullcline4(x)
    A4 = Complex((r/f)*((k*ϕ*x[2]/(F + kmin*x[1] - k*x[2]*r/(r + f*x[1]^2)) - 1)))^0.5 # only valid if ϕ = q/k
    erA4 = Inf
    if imag(A4) != 0
        return(erA4)
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
end

# jacobian should think a bit about how I do this
# x[1] = A, x[2] = B, x[3] = S, x[4] = W
function j!(J, x)
    J[1, 1] = -1
    J[1, 2] = -2*f*x[2]*k*x[3]*r/((kmin + K)*(r + f*x[2]^2)^2)
    J[1, 3] = k*r/((kmin + K)*(r + f*x[2]^2))
    J[1, 4] = Kmin/(kmin + K)
    J[2, 1] = -1
    denom = (kmin + K)*x[2] - Kmin*x[4]
    A = -((r/f)^0.5)*k*x[3]/(2*(denom)*((k*x[3]/(denom) - 1)^2))
    J[2, 2] = -(kmin + K)/(denom)
    J[2, 3] = A
    J[2, 4] = Kmin*A/(denom)
    J[3, 1] = -1
    J[3, 2] = -ϕ
    J[3, 3] = 0
    J[3, 4] = (Kmin/K)*(1 + ϕ)
    J[4, 1] = -1
    denom2 = F + kmin*x[2] - r/(r + f*x[2]^2)
    B = ((r/f)^0.5)*ϕ*k/(2*denom2*((ϕ*k*x[3]/(denom2) - 1)^0.5))
    J[4, 2] = B*x[3]*(kmin + 2*r*f*x[2]/((r + f*x[2]^2)^2)/denom2
    J[4, 3] = B
    J[4, 4] = 0
end



function find_zeros()
    # system of 6 equations to find the minimum of
    xs12 = nlsolve(f!, j!, [0; 0; 0; 0])
    print(xs12)
    print("\n")

end

@time find_zeros()
