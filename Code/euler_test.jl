#!/usr/bin/env julia
# euler_test.jl
# A script to run a euler integration of the equations with noise, adapted heavily from a
# matlab script Robert wrote
using DifferentialEquations
using Plots
import GR

# numerically calculates toggle switch
const Ω = 10 # volume
const k = 100
const K = k/Ω
const q = 10 #  asymmetric switches - second switch slower but same steady states.
const Q = q/Ω

const f = 0.01 # binding/unbinding
const r = 0.01
const a = 2 # Hill coefficient
const b = 2

function f1(du, u, p, t)
    du[1] = k*r/(r + f*u[2]^a) - K*u[1]
    du[2] = q*r/(r + f*u[1]^b) - Q*u[2]
end

function g1(du, u, p, t)
    du[1] = sqrt(abs((q*r/(r + f*u[1]^b) + Q*u[2]))/Ω)
    du[2] = sqrt(abs((k*r/(r + f*u[2]^a) + K*u[1]))/Ω)
end

function main()
    u0 = [10.0; 10.0]
    dt = (1/2)^(10)
    tspan = (0.0,50.0)
    prob = SDEProblem(f1, g1, u0, tspan) # SDEProblem
    sol = solve(prob, EM(), dt = dt)
    gr()
    plot(sol)
    savefig("../Results/Graph.png")
end

# run the main function
@time main()
