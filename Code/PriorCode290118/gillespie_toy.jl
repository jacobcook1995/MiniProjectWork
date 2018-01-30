#!/usr/bin/env julia
# Julia Script for Gillespie simulation of toy driven ecosystem

#####################################################################
# Main changes I've made to script are that it no longer cats each loop but
# instead cats every specified number of loops given by batchsize, this speeds
# the code up significantly. The production of the histogram is now done by 
# using the vectors to calculate the time spent at each level and then plotting
# these cumlative times as bar graphs. Otherwise the only change is the plotting 
# but this is mainly due to the different plotting syntax of julia
####################################################################

# Add package for plotting
using Plots
    
# Parameters
k = 10000 # steady state for A=k/K=1
K = 100 # K=k'
q = 100 # steady state for B=q/Q=1
Q = 1 # Q=q'
l = 0.001 # switching
batchsize = 1000000
 
Ti = 0  # initial time
Tf = 2000 # end of simulation time in s

At = []
Bt = []
T = []
 
Atemp = fill(0.0, batchsize) # preallocate
Btemp = fill(0.0, batchsize)
Ttemp = fill(0.0, batchsize)

Atemp[1] = 1
Btemp[1] = 1
Ttemp[1] = Ti
 
t = Ttemp[1] # initialise
A = Atemp[1]
B = Btemp[1]
i = 1

# Main loop
while t <= Tf
    i += 1
    # rates
    rates = [k, K*A, q, Q*B, l*A, l*B]
    r = rates/sum(rates)
    
    rone = rand() #first random number
   
    # which reaction?
    if rone < r[1]
        A += 1
    elseif rone < (r[1]+r[2])
        A += -1
    elseif rone < (r[1]+r[2]+r[3])
        B += 1
    elseif rone < (r[1]+r[2]+r[3]+r[4])
        B += -1
    elseif rone < (r[1]+r[2]+r[3]+r[4]+r[5])
        A += -1; B += 1
    else
        B += -1; A += 1
    end
    
    # update numbers
    Atemp[i] = A
    Btemp[i] = B
    
    # update time
    rtwo = rand()  # second random number
    t = t - log(rtwo)/sum(rates)
    Ttemp[i] = t
    if i == batchsize
        At = vcat(At, Atemp) # these three lines are probably the biggest use of time in this program
        Bt = vcat(Bt, Btemp)
        T = vcat(T, Ttemp)
        i = 0
        Atemp = fill(0.0, batchsize) # repreallocate
        Btemp = fill(0.0, batchsize)
        Ttemp = fill(0.0, batchsize)
    end
    
end

# Need to now cat the final incomplete batch, will test before I do this though
Averytemp = fill(0.0, i)
Bverytemp = fill(0.0, i)
Tverytemp = fill(0.0, i)
# Select the valued entries
for j = 1:i
    Averytemp[j] = Atemp[j]
    Bverytemp[j] = Btemp[j]
    Tverytemp[j] = Ttemp[j]
end
# Then cat onto vector
At = vcat(At, Averytemp)
Bt = vcat(Bt, Bverytemp)
T = vcat(T, Tverytemp)

# now rescale At and Bt based on the time interval they are valid for
maxA = convert(Int64,maximum(At)+1)
Atime = fill(0.0, maxA)
for j = 1:maxA
    for k = 1:length(At)
        if (At[k] == j-1)
            if k != length(At)
                Atime[j] += T[k+1] - T[k]
            else
                Atime[j] += Tf - T[k]
            end
        end
    end
end
maxB = convert(Int64,maximum(Bt)+1)
Btime = fill(0.0, maxB)
for j = 1:maxB
    for k = 1:length(Bt)
        if (Bt[k] == j-1)
            if k != length(Bt)
                Btime[j] += T[k+1] - T[k]
            else
                Btime[j] += Tf - T[k]
            end
        end
    end
end
# Finally calculate the moments
momentA = 0
for j = 1:length(Atime)
    momentA += (j-1)*Atime[j]/(Tf-Ti)
end
momentB = 0
for j = 1:length(Btime)
    momentB += (j-1)*Btime[j]/(Tf-Ti)
end

print("Gillespie done!\n")
 
gr() # set plotting back end to plotlyjs()
# Plot
pone = plot(T, [At Bt], xlabel = "Time", ylabel = "Solution", legend = false)
print("First plot done!\n")

# 1st histograms
binsA = 0:length(Atime)-1
ptwo = plot(binsA, Atime, xlabel = "A", ylabel = "Frequency", color = :blue, linetype=:bar, xlim = (-1,300), legend = false)
annotate!(150,20,text("<A>=$(momentA)",:left))
print("Second plot done!\n")
# 2nd histogram
binsB = 0:length(Btime)-1
pthree = plot(binsB, Btime, xlabel = "B", ylabel = "Frequency", color = :red, linetype=:bar, xlim = (-1,300), legend = false)
annotate!(150,20,text("<B>=$(momentB)",:left))
print("Third plot done!\n")

plot(pone, ptwo, pthree, layout=(3,1))#pthree

print("Plots combined!\n")

savefig("../Results/SolutionvsTime.png")
print("Plot Saved!\n")

# print moments to screen
print(momentA)
print("\n")
print(momentB)
print("\n")

print("Finished!\n")
