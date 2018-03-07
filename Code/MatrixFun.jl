#!/usr/bin/env julia
# MatrixFun.jl
# A script to test matrix algebra in julia hopefully developing to the point
# where I can write a function to calculate the diffusion matrix for arbitary noise

using Plots
using Roots
using SymPy
import GR # Need this to stop world age plotting error?

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
const F = 250 # removal rate

# Multiplicative Guassian noise matrix
# A = x[1], B = x[2], W = x[3], S = x[4]
function e!(E, x)
    E[1,2:4] = E[2,1] = E[2,3:4] = E[3,4] = E[4,3] = 0
    E[1,1] = sqrt(k*x[4]*r/(r + f*x[2]^2) + K*x[1] + kmin*x[1] + Kmin*x[3]) #gA
    E[2,2] = sqrt(q*x[4]*r/(r + f*x[1]^2) + Q*x[2] + qmin*x[2] + Qmin*x[3]) #gB
    E[3,1] = -sqrt(K*x[1] + Kmin*x[3]) #-gWA
    E[3,2] = -sqrt(Q*x[2] + Qmin*x[3]) #-gWB
    E[3,3] = sqrt(F) #gW
    E[4,1] = -sqrt(k*x[4]*r/(r + f*x[2]^2) + kmin*x[1]) #-gSA
    E[4,2] = -sqrt(q*x[4]*r/(r + f*x[1]^2) + qmin*x[2]) #-gSB
    E[4,4] = sqrt(F) #gS
    return E
end

# Gonna try numerically first
function numerically()
    e = [ 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 ]
    e = e!(e, [ 10; 10; 10; 10 ])
    eT = transpose(e)
    D = e*eT
    H = D*D
    Dmin = inv(H)
    return(Dmin)
end

# Now gonna try symbolically
function symbolically()
    # # Locally overwrites global constants to be variables
    r, f, q, qmin, Q, Qmin, k, kmin, K, Kmin, F = symbols("r,f,q,qmin,Q,Qmin,k,kmin,K,Kmin,F")
    A, B, W, S = symbols("A,B,W,S")
    #
    # # Make a symbolic version of the matrix, needs no input in this case
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

    # Dmin = inv(D)
    K1 = 10
    k1 = K1*Ω # steady state for A=k/K=1
    Q1 = K1*ϕ
    q1 = Q1*Ω
    kmin1 = 10.0^-20 # set all too 10.0^-20 for now
    Kmin1 = 10.0^-20
    qmin1 = 10.0^-20
    Qmin1 = 10.0^-20
    f1 = 1000/(Ω^2) # Promoter switching
    r1 = 10
    F1 = 250 # removal rate
    eT = transpose(e)
    D = e*eT
    #Dmin = inv(D)
    H = D*D
    Dmin = inv(H)
    Dmin = subs(Dmin, A, 10)
    Dmin = subs(Dmin, B, 10)
    Dmin = subs(Dmin, S, 10)
    Dmin = subs(Dmin, W, 10)
    Dmin = subs(Dmin, K, K1)
    Dmin = subs(Dmin, k, k1)
    Dmin = subs(Dmin, Q, Q1)
    Dmin = subs(Dmin, q, q1)
    Dmin = subs(Dmin, kmin, kmin1)
    Dmin = subs(Dmin, Kmin, Kmin1)
    Dmin = subs(Dmin, qmin, qmin1)
    Dmin = subs(Dmin, Qmin, Qmin1)
    Dmin = subs(Dmin, f, f1)
    Dmin = subs(Dmin, r, r1)
    Dmin = subs(Dmin, F, F1) |> float
    return(Dmin)
end

@time Dmin1 = numerically()
@time Dmin2 = symbolically()
#
print(Dmin1)
print("\n")
print(Dmin2)
print("\n")
