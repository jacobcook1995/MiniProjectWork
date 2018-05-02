#!/usr/bin/env julia
# euler_test.jl
# A script to run a euler integration of the equations with noise, adapted heavily from a
# matlab script Robert wrote
using DifferentialEquations
using Plots
using Roots
using SymPy
import GR

# Parameters
const Ω = 300 # system size
const ϕ = 0.1 # ratio ϕ = q/k
const K = 10
const k = K*Ω # steady state for A=k/K=1
const Q = K*ϕ
const q = Q*Ω
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const Kmin = 10.0^-20
const qmin = 10.0^-20
const Qmin = 10.0^-20
const f = 1000/(Ω^2) # Promoter switching
const r = 10
const F = 300 # removal rate
const Ne = 12000 # number of elements in the system
const high2low = false

# Diffusion matrix function in inverse form, this will become a global constant matrix
function Dmin1()
    A, B, W, S = symbols("A,B,W,S")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{Sym}(4,4)
    e[1,2:4] = e[2,1] = e[2,3:4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + K*A + kmin*A + Kmin*W) #gA
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + Q*B + qmin*B + Qmin*W) #gB
    e[3,1] = -sqrt(K*A + Kmin*W) #-gWA
    e[3,2] = -sqrt(Q*B + Qmin*W) #-gWB
    e[3,3] = sqrt(F) #gW
    e[4,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A) #-gSA
    e[4,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B) #-gSB
    e[4,4] = sqrt(F) #gS
    # Now do the transformations required
    eT = transpose(e)
    D = e*eT
    Dmin = inv(D)
    return(Dmin)
end

# Inverse Diffusion matrix containing the noise on each term (squared)
function D!(D, x)
    if x[1] < 0 || x[2] < 0 || x[3] < 0 || x[4] < 0
        print("$(x[1]),$(x[2]),$(x[3]),$(x[4])\n")
    end
    A, B, W, S = symbols("A,B,W,S")
    D = subs(Dminconst, A, x[1]) |> Sym
    D = subs(D, B, x[2]) |> Sym
    D = subs(D, W, x[4]) |> Sym
    D = subs(D, S, x[4]) |> float
    return(D)
end

# x[1] = A, x[2] = B, x[3] = W, x[4] = S
function f!(F1, x)
    F1[1] = k*x[4]*r/(r + f*x[2]*x[2]) - (K + kmin)*x[1] + Kmin*x[3]
    F1[2] = q*x[4]*r/(r + f*x[1]*x[1]) - (Q + qmin)*x[2] + Qmin*x[3]
    F1[3] = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3]
    F1[4] = -k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] - q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2]
    return(F1)
end

function fprim!(Fprim1, x)
    Fprim1[1] = -(K + kmin)
    Fprim1[2] = -(Q + qmin)
    Fprim1[3] = -Kmin*(1 + ϕ)
    Fprim1[4] = -k*r/(r + f*x[2]*x[2]) - q*r/(r + f*x[1]*x[1])
    return(Fprim1)
end

function Fσ!(Fσ)
    Fσ[1] = 0
    Fσ[2] = 0
    Fσ[3] = F
    Fσ[4] = -F
    return(Fσ)
end

# funtion to find the zeros of the system should find 3 zeros ss1, ss2 and the saddle point
function Zeros()
    g(x) = F./(1 + ((r + f*x.^2)./(ϕ*r + (f/ϕ)*((F/K - x).^2)))) + K.*x - F
    three = false
    n = 0
    As = []
    while three == false
        As = fzeros(g, 0, 2*F/K, order = 1)
        n = length(As)
        if n == 3
            three = true
        end
    end
    Bs = zeros(n)
    Ss = zeros(n)
    Ws = zeros(n)
    for i = 1:n
        Bs[i] = (1/ϕ)*((F/K) - As[i])
        Ss[i] = F/(k*r*(1/(r + f*Bs[i]^2) + ϕ/(r + f*As[i]^2)))
        Ws[i] = Ne - As[i] - Bs[i] - Ss[i]
    end
    sad = [ As[2]; Bs[2]; Ws[2]; Ss[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ws[1]; Ss[1] ]
        ss2 = [ As[3]; Bs[3]; Ws[3]; Ss[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ws[3]; Ss[3] ]
        ss2 = [ As[1]; Bs[1]; Ws[1]; Ss[1] ]
    end
    print(ss1)
    print("\n")
    print(sad)
    print("\n")
    print(ss2)
    print("\n")
    return(ss1,sad,ss2)
end

function f1(du, u, p, t)
    du[1] = k*u[4]*r/(r + f*u[2]^2) - (K + kmin)*u[1] + Kmin*u[3]
    du[2] = q*u[4]*r/(r + f*u[1]^2) - (Q + qmin)*u[2] + Qmin*u[3]
    du[3] = K*(u[1] + ϕ*u[2]) - Kmin*(1 + ϕ)*u[3] - F
    du[4] = -k*u[4]*r/(r + f*u[2]^2) + kmin*u[1] - q*u[4]*r/(r + f*u[1]^2) + qmin*u[2] + F
end

function g1(du, u, p, t)
    # Need to define a 4 by four matrix as I have 4 weiner processes and 4 dependant variables
    du[1,2:4] = du[2,1] = du[2,3:4] = du[3,4] = du[4,3] = 0
    du[1,1] = sqrt((k*u[4]*r/(r + f*u[2]^2) + (K + kmin)*u[1] + Kmin*u[3])/Ω)
    du[2,2] = sqrt((q*u[4]*r/(r + f*u[1]^2) + (Q + qmin)*u[2] + Qmin*u[3])/Ω)
    du[3,3] = sqrt((F)/Ω)
    du[4,4] = sqrt((F)/Ω)
    du[3,1] = -sqrt((K*u[1] + Kmin*u[3])/Ω)
    du[3,2] = -sqrt((Q*u[2] + Qmin*u[3])/Ω)
    du[4,1] = -sqrt((k*u[4]*r/(r + f*u[2]^2) + kmin*u[1])/Ω)
    du[4,2] = -sqrt((q*u[4]*r/(r + f*u[1]^2) + qmin*u[2])/Ω)
end

function main()
    ss1, sad, ss2 = Zeros()
    # First generate a run of the problem
    u₀ = ss1
    dt = (1/2)^(12)
    tspan = (0.0,15)
    prob = SDEProblem(f1, g1, u₀, tspan, noise_rate_prototype = zeros(4,4)) # SDEProblem
    sol = DifferentialEquations.solve(prob, EM(), dt = dt) # To avoid namespace conflict with SymPy
    # Make temporary holders for variables at each step
    Dt = zeros(4,4)
    ft = zeros(4)
    thivt = zeros(4)
    Fσt = zeros(4)
    Fprimt = zeros(4)
    # Start as scalers can extend to arrays or vectors
    F1 = 0 # velocity^2 # kinetics
    F2 = 0 # force^2 # kinetics
    F3 = 0 # Constant^2 # kinetics
    F4 = 0 # entprod
    F5 = 0 # Velocity*constant # suppossedly zero
    F6 = 0 # constant*force # kinetics
    F7 = 0 # curvature
    # Same but cross terms
    F11 = 0 # velocity^2 # kinetics
    F12 = 0 # force^2 # kinetics
    F13 = 0 # Constant^2 # kinetics
    F14 = 0 # entprod
    F15 = 0 # Velocity*constant # suppossedly zero
    F16 = 0 # constant*force # kinetics

    # Now integrate over the last 50 steps to find the behaviour of all the different contributions
    for j = 1:60000
        Dt = D!(Dt, sol[:,end+1-j])
        ft = f!(ft, sol[:,end+1-j])
        thivt = (sol[:,end+1-j] - sol[:,end-j])/(dt)
        Fσt = Fσ!(Fσt)
        Fprimt = fprim!(Fprimt, sol[:,end+1-j])
        for l = 1:4
            for m = 1:4
                if l == m
                    F1 += 0.5*thivt[l]*Dt[l,m]*thivt[m]*dt
                    F2 += 0.5*ft[l]*Dt[l,m]*ft[m]*dt
                    F3 += 0.5*Fσt[l]*Dt[l,m]*Fσt[m]*dt
                    F4 -= thivt[l]*Dt[l,m]*ft[m]*dt
                    F5 -= thivt[l]*Dt[l,m]*Fσt[m]*dt
                    F6 += Fσt[l]*Dt[l,m]*ft[m]*dt
                else
                    F11 += 0.5*thivt[l]*Dt[l,m]*thivt[m]*dt
                    F12 += 0.5*ft[l]*Dt[l,m]*ft[m]*dt
                    F13 += 0.5*Fσt[l]*Dt[l,m]*Fσt[m]*dt
                    F14 -= thivt[l]*Dt[l,m]*ft[m]*dt
                    F15 -= thivt[l]*Dt[l,m]*Fσt[m]*dt
                    F16 += Fσt[l]*Dt[l,m]*ft[m]*dt
                end

            end
            F7 -= 0.5*Fprimt[l]*dt/Ω
        end
    end
    print("$(F1+F11)\n")
    print("$(F2+F12)\n")
    print("$(F3+F13)\n")
    print("$(F4+F14)\n")
    print("$(F5+F15)\n")
    print("$(F6+F16)\n")
    print("$(F7)\n")
    print("Kinetics (Proper) = $(F1+F2+F3+F6)\n")
    print("Entropy Production (Proper) = $(F4)\n")
    print("Zero (Proper) = $(F5)\n")
    print("Curvature = $(F7)\n")
    print("Kinetics (Cross) = $(F11+F12+F13+F16)\n")
    print("Entropy Production (Cross) = $(F14)\n")
    print("Zero (Cross) = $(F15)\n")
end


# Symbolic matrix set as global constant
D = Dmin1()
const Dminconst = D

# run the main function
# @time main()
# @time test()
