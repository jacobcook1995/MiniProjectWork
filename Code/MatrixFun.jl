#!/usr/bin/env julia
# MatrixFun.jl
# A script to test matrix algebra in julia hopefully developing to the point
# where I can write a function to calculate the diffusion matrix for arbitary noise

using Plots
using Roots
using SymPy
import GR # Need this to stop world age plotting error?

# Parameters
const Ω = 30 # system size
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
function e(x)
    E = zeros(4,4)
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
    h = 1*(10.0^-1.5)
    diff = h*[ 1; 0; 0; 0 ]
    v = [ 10; 10; 10; 10 ]
    v1 = v + diff
    v2 = v + 2*diff
    vmin1 = v - diff
    vmin2 = v - 2*diff
    # use higher order numeric differentiation here
    e1 = e(v1)
    eT1 = transpose(e1)
    D1 = e1*eT1
    e2 = e(v2)
    eT2 = transpose(e2)
    D2 = e2*eT2
    emin1 = e(vmin1)
    eTmin1 = transpose(emin1)
    Dmin1 = emin1*eTmin1
    emin2 = e(vmin2)
    eTmin2 = transpose(emin2)
    Dmin2 = emin2*eTmin2
    Dinv1 = inv(D1)
    Dinv2 = inv(D2)
    Dinvmin1 = inv(Dmin1)
    Dinvmin2 = inv(Dmin2)
    Dminprim = (-Dinv2 + 8*Dinv1 - 8*Dinvmin1 + Dinvmin2)/(12*h)
    return(Dminprim)
end

# Now gonna try symbolically
function symbolically1()
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
    Dmin1 = inv(D)
    Dmin1prim = diff(Dmin1, A)
    Dmin1prim = subs(Dmin1prim, K, K1) |> Sym
    Dmin1prim = subs(Dmin1prim, k, k1) |> Sym
    Dmin1prim = subs(Dmin1prim, Q, Q1) |> Sym
    Dmin1prim = subs(Dmin1prim, q, q1) |> Sym
    Dmin1prim = subs(Dmin1prim, kmin, kmin1) |> Sym
    Dmin1prim = subs(Dmin1prim, Kmin, Kmin1) |> Sym
    Dmin1prim = subs(Dmin1prim, qmin, qmin1) |> Sym
    Dmin1prim = subs(Dmin1prim, Qmin, Qmin1) |> Sym
    Dmin1prim = subs(Dmin1prim, f, f1) |> Sym
    Dmin1prim = subs(Dmin1prim, r, r1) |> Sym
    Dmin1prim = subs(Dmin1prim, F, F1) |> Sym
    return(Dmin1prim)
end

# Now gonna try symbolically
function symbolically2()
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

    eT = transpose(e)
    D = e*eT
    Dmin1 = inv(D)
    Dmin1prim = diff(Dmin1, A)
    return(Dmin1prim)
end

D1 = symbolically1()
const D1c = D1
D2 = symbolically2()
const D2c = D2
function sub1()
    A, B, W, S = symbols("A,B,W,S")
    D3 = subs(D1c, A, 10.7) |> Sym
    D3 = subs(D3, B, 103) |> Sym
    D3 = subs(D3, W, 10000) |> Sym
    D3 = subs(D3, S, 26) |> float
    return(D3)
end

function sub2()
    A, B, W, S = symbols("A,B,W,S")
    r, f, q, qmin, Q, Qmin, k, kmin, K, Kmin, F = symbols("r,f,q,qmin,Q,Qmin,k,kmin,K,Kmin,F")
    # Now do the transformations required
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
    D4 = subs(D2c, K, K1) |> Sym
    D4 = subs(D4, A, 10.7) |> Sym
    D4 = subs(D4, B, 103) |> Sym
    D4 = subs(D4, F, F1) |> Sym
    D4 = subs(D4, Kmin, Kmin1) |> Sym
    D4 = subs(D4, Q, Q1) |> Sym
    D4 = subs(D4, Qmin, Qmin1) |> Sym
    D4 = subs(D4, S, 26) |> Sym
    D4 = subs(D4, W, 10000) |> Sym
    D4 = subs(D4, f, f1) |> Sym
    D4 = subs(D4, k, k1) |> Sym
    D4 = subs(D4, kmin, kmin1) |> Sym
    D4 = subs(D4, q, q1) |> Sym
    D4 = subs(D4, qmin, qmin1) |> Sym
    D4 = subs(D4, r, r1) |> float
    return(D4)
end


D4 = sub2()
D3 = sub1()
print("$D3\n")
print("$D4\n")
print(D3-D4)
print("\n")
