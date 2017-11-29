#!/usr/bin/env julia
# Julia Script for Gillespie simulation of toy driven ecosystem

# Add package for plotting
using Plots
    
# Parameters
k = 100 # steady state for A=k/K=1
K = 100 # K=k'
q = 1 # steady state for B=q/Q=1
Q = 1 # Q=q'
l = 0.001 # switching
 
Ti = 0  # initial time
Tf = 1000 # end of simulation time in s
 
At=[1] # trajectories
Bt=[1]
T=[Ti]
 
t=T[1] # initialise
A=At[1]
B=Bt[1]

# Main loop
while t <= Tf
    # rates
    rates=[k, K*A, q, Q*B, l*A, l*B]
    r=rates/sum(rates)
    
    rone=rand(); #first random number
   
    # which reaction?
    if rone<r[1]
        A = A + 1
    elseif r[1]<=rone && rone<r[1]+r[2]
        A = A - 1
    elseif r[1]+r[2]<=rone && rone<r[1]+r[2]+r[3]
        B = B + 1
    elseif r[1]+r[2]+r[3]<=rone && rone<r[1]+r[2]+r[3]+r[4]
        B = B - 1
    elseif r[1]+r[2]+r[3]+r[4]<=rone && rone<r[1]+r[2]+r[3]+r[4]+r[5]
        A = A - 1; B = B + 1
    else
        B = B - 1; A = A + 1
    end
    
        
    # update numbers
    At = vcat(At,A)
    Bt = vcat(Bt,B)
    
    # update time
    rtwo=rand()  # second random number
    t = t - log(rtwo)/sum(rates)
    T = vcat(T,t)
    
end
 
plotlyjs() # set plotting back end to plotlyjs()
# Plot
pone = plot(T, At, label = "A(t)", xlabel = "Time", ylabel = "Solution")
pone = plot!(pone, T, Bt, label = "B(t)") # mutate plot

# 1st histograms
ptwo = histogram(At, bins = 0:maximum(At), xlabel = "A", ylabel = "Frequency", color = :blue)
# 2nd histogram
pthree = histogram(Bt, bins = 0:maximum(Bt), xlabel = "B", ylabel = "Frequency", color = :red)
plot(pone,ptwo,pthree,layout=(3,1))

savefig( "../Results/SolutionvsTime.png")

# Now calculate the moments
momentA = sum(At)/length(At)
momentB = sum(Bt)/length(Bt)
print(momentA)
print("\n")
print(momentB)
print("\n")

print("Finished!\n")
