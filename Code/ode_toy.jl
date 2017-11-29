#!/usr/bin/env julia
# Julia Script for ODE model of toy driven ecosystem

# First the includes
using DifferentialEquations
using Plots

# Then function definitions
function toy(t,u,du)
    # Constants
    k=100 # steady state for A=k/K=1
    K=100 # K=k'
    q=1 # steady state for B=q/Q=1
    Q=1 # Q=q'
    l=0.001 #switching
    du[1] = k + l*u[2] - (l+K)*u[1]
    du[2] = q + l*u[1] - (l+Q)*u[2]
end


# Then body of code
# 1. time course
prob = ODEProblem(toy,[0.0;1.0],(0.0,1.0))
sol = solve(prob) # started on this, not clear if I'm right or not

plotlyjs() # set plotting back end to plotlyjs()
# Plot
p = plot(sol, vars=(0,1), label = "A(t)", xlabel = "Time", ylabel = "#", color = :red)
p = plot!(p,sol, vars=(0,2), label = "B(t)", color = :blue, title = "Solution of driven bistable system")
savefig( "../Results/SolutionvsTimeODE.png")
 
# 2. steady state values from analytical theory

# Constants
k=100 # steady state for A=k/K=1
K=100 # K=k'
q=1 # steady state for B=q/Q=1
Q=1 # Q=q'
l=0.001 #switching 
A = (q*l+(l+Q)*k)/((l+Q)*(l+K)-l*l)
B = (l+K)/l*(q*l+(l+Q)*k)/((l+Q)*(l+K)-l*l)-k/l
print(A)
print("\n")
print(B)
print("\n")

print("Finished\n")
