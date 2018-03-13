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
const high2low = true

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
    A, B, W, S = symbols("A,B,W,S")
    D = subs(Dminconst, A, x[1]) |> Sym
    D = subs(D, B, x[2]) |> Sym
    D = subs(D, W, x[4]) |> Sym
    D = subs(D, S, x[4]) |> float
    return D
end

# x[1] = A, x[2] = B, x[3] = W, x[4] = S
function f!(F1, x)
    F1[1] = k*x[4]*r/(r + f*x[2]*x[2]) - (K + kmin)*x[1] + Kmin*x[3]
    F1[2] = q*x[4]*r/(r + f*x[1]*x[1]) - (Q + qmin)*x[2] + Qmin*x[3]
    F1[3] = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3] - F
    F1[4] = -k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] - q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F
    return F1
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
    return (ss1,sad,ss2)
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
    du[1,1] = sqrt(abs(k*u[4]*r/(r + f*u[2]^2) + (K + kmin)*u[1] + Kmin*u[3])/Ω)
    du[2,2] = sqrt(abs(q*u[4]*r/(r + f*u[1]^2) + (Q + qmin)*u[2] + Qmin*u[3])/Ω)
    du[3,3] = sqrt(abs(F)/Ω)
    du[4,4] = sqrt(abs(F)/Ω)
    du[3,1] = -sqrt(abs(K*u[1] + Kmin*u[3])/Ω)
    du[3,2] = -sqrt(abs(Q*u[2] + Qmin*u[3])/Ω)
    du[4,1] = -sqrt(abs(k*u[4]*r/(r + f*u[2]^2) + kmin*u[1])/Ω)
    du[4,2] = -sqrt(abs(q*u[4]*r/(r + f*u[1]^2) + qmin*u[2])/Ω)
end

function main()
    ss1, sad, ss2 = Zeros()
    Dt = zeros(4,4)
    ht = zeros(4)
    D = zeros(4,4)
    h = zeros(4)
    thiv = zeros(4)
    states = [ ss1 ]

    # first component gives states, third component gives i, fourth gives j,
    A = zeros(2,4,4,4) # second gives component [KE, PE, EP1, EP2]
    for i = 1:length(states)
        u₀ = states[i]
        dt = (1/2)^(10)
        tspan = (0.0,50.0)
        prob = SDEProblem(f1, g1, u₀, tspan, noise_rate_prototype = zeros(4,4)) # SDEProblem
        sol = DifferentialEquations.solve(prob, EM(), dt = dt) # To avoid namespace conflict with SymPy
        for j = 1:5000
            Dt = D!(Dt, sol[:,end+1-j])
            ht = f!(ht, sol[:,end+1-j])
            thivt = (sol[:,end+1-j] - sol[:,end-j])/(dt)
            D += Dt
            h += ht
            thiv += thivt
        end
        plot([sol[1,:],sol[2,:]])
        savefig("../Results/Graph$(i).png")
        for j = 1:4
            for l = 1:4
                for m = 1:4
                    if j == 1
                        A[i,j,l,m] = thiv[l]*D[l,m]*thiv[m]
                    elseif j == 2
                        A[i,j,l,m] = h[l]*D[l,m]*h[m]
                    elseif j == 3
                        A[i,j,l,m] = -h[l]*D[l,m]*thiv[m]
                    else
                        A[i,j,l,m] = -thiv[l]*D[l,m]*h[m]
                    end
                end
            end
        end
    end
print(A[1,1,:,:])
print("\n")
print(A[1,2,:,:])
print("\n")
print(A[1,3,:,:])
print("\n")
print(A[1,4,:,:])
print("\n")

end

# Symbolic matrix set as global constant
D = Dmin1()
const Dminconst = D

# run the main function
@time main()
