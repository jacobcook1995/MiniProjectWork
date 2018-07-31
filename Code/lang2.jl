#!/usr/bin/env julia
# lang2.jl
# A script to efficently perform a langevin calculation of the reduced 2 species model
# This script then tries to calculate the entropy production and action from
# the generated Langevin trajectories

using DifferentialEquations
using Plots
using Roots
import GR

# finction to find start end and saddle points
function nullcline(r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64)
    a = 2
    b = 2
    A1(x) = k*r/(K*(r+f*x^a))
    A2(x) = (r/f*(q/(Q*x)-1))^(1/b)
    g(x) = k*r/(K*(r+f*x^a)) - (r/f*(q/(Q*x)-1))^(1/b)
    xs = fzeros(g, 0, q/Q)
    ss1 = [A1(xs[1]); xs[1]]
    sad = [A1(xs[2]); xs[2]]
    ss2 = [A1(xs[3]); xs[3]]
    println(ss1)
    println(sad)
    println(ss2)
    return (ss1,sad,ss2)
end

# function describing bistable system
# p = [ Ω, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r]
#       1, 2, 3, 4, 5,    6,    7,    8,    9, 10,11]
function bis(du,u,p,t)
    du[1] = p[3]*p[11]/(p[11] + p[10]*(u[2]^2)) - p[6]*u[1] - p[2]*u[1] + p[7]
    du[2] = p[5]*p[11]/(p[11] + p[10]*(u[1]^2)) - p[8]*u[2] - p[4]*u[2] + p[9]
    return(du)
end
# diagonal noise function
function σ_bis(du,u,p,t)
    du[1] = (p[3]*p[11]/(p[11] + p[10]*(u[2]^2)) + p[6]*u[1] + p[2]*u[1] + p[7])/sqrt(p[1])
    du[2] = (p[5]*p[11]/(p[11] + p[10]*(u[1]^2)) + p[8]*u[2] + p[4]*u[2] + p[9])/sqrt(p[1])
    return(du)
end
# function to return fs
function fs(u,p)
    fA = p[3]*p[11]/(p[11] + p[10]*(u[2]^2)) - p[6]*u[1] - p[2]*u[1] + p[7]
    fB = p[5]*p[11]/(p[11] + p[10]*(u[1]^2)) - p[8]*u[2] - p[4]*u[2] + p[9]
    return(fA,fB)
end
function fsp(u,p)
    fpA = -p[6] - p[2]
    fpB = -p[8] - p[4]
    return(fpA,fpB)
end
# function to return Ds
function Ds(u,p)
    DA = p[3]*p[11]/(p[11] + p[10]*(u[2]^2)) + p[6]*u[1] + p[2]*u[1] + p[7]
    DB = p[5]*p[11]/(p[11] + p[10]*(u[1]^2)) + p[8]*u[2] + p[4]*u[2] + p[9]
    return(DA,DB)
end
# function to generate path entropy poroduction
function Entp(pos::Array{Array{Float64,1},1},times::Array{Float64,1},p::Array{Float64,1})
    n = length(times) - 1
    A = zeros(n)
    qq = zeros(n)
    qf = zeros(n)
    ff = zeros(n)
    for i = 1:n
        τ = times[i+1] - times[i]
        posA = pos[i][1]#(pos[i+1][1] + pos[i][1])/2#
        posB = pos[i][2]#(pos[i+1][2] + pos[i][2])/2#
        fA, fB = fs([posA, posB], p)
        DA, DB = Ds([posA, posB], p)
        qdA = (pos[i+1][1] - pos[i][1])/τ
        qdB = (pos[i+1][2] - pos[i][2])/τ
        qq[i] = ((qdA^2)/DA + (qdB^2)/DB)*τ
        qf[i] = ((fA*qdA)/DA + (fB*qdB)/DB)*τ
        ff[i] = ((fA^2)/DA + (fB^2)/DB)*τ
        A[i] = 0.5*(qq[i]+ff[i]) - qf[i]
    end
    return(A,qq,qf,ff)
end

function main()
    # General parameters
    Ω = 150.0
    K = 10.0
    k = K*Ω # steady state for A=k/K=1
    Q = 1.0
    q = Q*Ω
    kmin = 10.0^-20 # set all too 10.0^-20 for now
    Kmin = (10.0^-20)*Ω
    qmin = 10.0^-20
    Qmin = (10.0^-20)*Ω
    f = 1000.0/(Ω^2) # Promoter switching
    r = 10.0
    p = [ Ω, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r]
    star, mid, fin = nullcline(r,f,K,Q,k,q,kmin,qmin)
    prob_sde_bis = SDEProblem(bis,σ_bis,fin,(0.0,10.0),p)
    sol = solve(prob_sde_bis)
    # now need a function to calculate action and entropy prdouction along path
    tf = sol.t[end]
    println(length(sol.t))
    A, qq, qf, ff = Entp(sol.u,sol.t,p)
    FA, FB = fsp([0.0, 0.0], p)
    println(sum(qq)/tf)
    println(sum(qf)/tf)
    println(sum(ff)/tf)
    println(sum(A)*Ω/tf)
    println(sum(qf)*Ω/tf)
    println(FA+FB)
    plot(sol.t[1:end-1],qf)
    savefig("../Results/graph1.png")
    plot(sol,vars=(0,1))
    savefig("../Results/graph2.png")
    plot(sol,vars=(0,2))
    savefig("../Results/graph3.png")
    return(nothing)
end

@time main()
