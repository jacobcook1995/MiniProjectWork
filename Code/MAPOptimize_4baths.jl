#!/usr/bin/env julia
# MAPOptimize_4baths.jl
# A script to optimize the path function of the minimum action path in order
# to determine the minimum path action for the driven bistable gene switch with
# cooperativity.
# This time it is done for the case of 4 baths, A, B, S, W
# This work draws heavily from the work of Ruben Perez-Carrasco et al (2016)

# Putting the relevant imports in
using Optim
using Plots
using Roots
import GR # Need this to stop world age plotting error?

# Parameters
const Ω = 30 # system size
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
const Ne = 2000 # number of elements in the system

# Then set parameters of the optimization
const N = 150 # number of segments optimised over

# Now construct the three relevant vectors of equations
# x[1] = A, x[2] = B, x[3] = W, x[4] = S
function f!(F1, x)
    F1[1] = k*x[4]*r/(r + f*x[2]*x[2]) - (K + kmin)*x[1] + Kmin*x[3]
    F1[2] = q*x[4]*r/(r + f*x[1]*x[1]) - (Q + qmin)*x[2] + Qmin*x[3]
    F1[3] = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3] - F
    F1[4] = -k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] - q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F
    return F1
end

function fA!(FA, x)
    FA[1] = -(K + kmin)
    FA[2] = -2*x[4]*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    FA[3] = K
    FA[4] = kmin + 2*x[4]*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    return FA
end

function fB!(FB, x)
    FB[1] = -2*x[4]*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    FB[2] = -(Q + qmin)
    FB[3] = Q
    FB[4] = qmin + 2*x[4]*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    return FB
end

function fW!(FW, x)
    FW[1] = Kmin
    FW[2] = Qmin
    FW[3] = -(1 + ϕ)*Kmin
    FW[4] = 0
    return FW
end

function fS!(FS, x)
    FS[1] = k*r/(r + f*x[2]*x[2])
    FS[2] = q*r/(r + f*x[1]*x[1])
    FS[3] = 0
    FS[4] = -k*r/(r + f*x[2]*x[2]) - q*r/(r + f*x[1]*x[1])
    return FS
end

# Diffusion matrix containing the noise on each term (squared)
function D!(D, x)
    D[1,1] = k*x[4]*r/(r + f*x[2]*x[2]) + (K + kmin)*x[1] + Kmin*x[3]
    D[2,2] = q*x[4]*r/(r + f*x[1]*x[1]) + (Q + qmin)*x[2] + Qmin*x[3]
    D[3,3] = K*(x[1] + ϕ*x[2]) + Kmin*(1 + ϕ)*x[3] + F
    D[4,4] = k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] + q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F
    D[1,2:4] = D[2,1] = D[2,3:4] = D[3,1:2] = D[3,4] = D[4,1:3] = 0
    return D
end

function D2!(D2, x)
    D2[1,1] = (k*x[4]*r/(r + f*x[2]*x[2]) + (K + kmin)*x[1] + Kmin*x[3])^2
    D2[2,2] = (q*x[4]*r/(r + f*x[1]*x[1]) + (Q + qmin)*x[2] + Qmin*x[3])^2
    D2[3,3] = (K*(x[1] + ϕ*x[2]) + Kmin*(1 + ϕ)*x[3] + F)^2
    D2[4,4] = (k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] + q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F)^2
    D2[1,2:4] = D2[2,1] = D2[2,3:4] = D2[3,1:2] = D2[3,4] = D2[4,1:3] = 0
    return D2
end

function DA!(DA, x)
    DA[1,1] = (K + kmin)
    DA[2,2] = -2*q*x[4]*r*x[1]/((r + f*x[1]*x[1])^2)
    DA[3,3] = K
    DA[4,4] = kmin - 2*q*x[4]*r*x[1]/((r + f*x[1]*x[1])^2)
    DA[1,2:4] = DA[2,1] = DA[2,3:4] = DA[3,1:2] = DA[3,4] = DA[4,1:3] = 0
    return DA
end

function DB!(DB, x)
    DB[1,1] = -2*k*x[4]*r*x[2]/((r + f*x[2]+x[2])^2)
    DB[2,2] = (Q + qmin)
    DB[3,3] = Q
    DB[4,4] = qmin - 2*k*x[4]*r*x[2]/((r + f*x[2]*x[2])^2)
    DB[1,2:4] = DB[2,1] = DB[2,3:4] = DB[3,1:2] = DB[3,4] = DB[4,1:3] = 0
    return DB
end

function DW!(DW, x)
    DW[1,1] = Kmin
    DW[2,2] = Qmin
    DW[3,3] = Kmin*(1 + ϕ)
    DW[4,4] = 0
    DW[1,2:4] = DW[2,1] = DW[2,3:4] = DW[3,1:2] = DW[3,4] = DW[4,1:3] = 0
    return DW
end

function DS!(DS, x)
    DS[1,1] = k*r/(r + f*x[2]*x[2])
    DS[2,2] = q*r/(r + f*x[1]*x[1])
    DS[3,3] = 0
    DS[4,4] = k*r/(r + f*x[2]*x[2]) + q*r/(r + f*x[1]*x[1])
    DS[1,2:4] = DS[2,1] = DS[2,3:4] = DS[3,1:2] = DS[3,4] = DS[4,1:3] = 0
    return DS
end
