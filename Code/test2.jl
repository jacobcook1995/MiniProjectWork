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
Ω = 30 # system size
k = 100 # steady state for A=k/K=1
K = k/Ω # K=k'
q = 10 # steady state for B=q/Q=1
Q = q/Ω # Q=q'
f = 1 # Promoter switching
rr = 1
batchsizeA = 1000000
batchsizeB = 1000000

Na = 0 # counts of A and B, mainly to give an idea of statisical weight later on
Nb = 0

Ti = 0.0  # initial time
Tf = 100000 # end of simulation time in s

At = []
Bt = []
TA = []
TB = []

Atemp = fill(0, batchsizeA) # preallocate
Btemp = fill(0, batchsizeB)
TtempA = fill(0.0, batchsizeA)
TtempB = fill(0.0, batchsizeB)

Atemp[1] = 1
Btemp[1] = Ω
TtempA[1] = Ti
TtempB[1] = Ti

t = TtempA[1] # initialise
A = Atemp[1]
B = Btemp[1]
i = 2 # start counts at two as the first element to fill is two
j = 2

a = 1
b = 1

# Main loop
while t <= Tf
    # rates
    rates = [K*A, a*k, Q*B, b*q, f*B*(B-1)*a, (1-a)*rr, f*A*(A-1)*b, (1-b)*rr]
    r = rates/sum(rates)

    rone = rand() #first random number

    updateA = 0
    updateB = 0

    # which reaction?
    if rone < r[1]
        A -= 1
        updateA = 1
    elseif rone < (r[1]+r[2])
        A += 1
        updateA = 1
    elseif rone < (r[1]+r[2]+r[3])
        B -= 1
        updateB = 1
    elseif rone < (r[1]+r[2]+r[3]+r[4])
        B += 1
        updateB = 1
    elseif rone < (r[1]+r[2]+r[3]+r[4]+r[5])
        a -= 1
        B -= 2 # With cooperativity
        updateB = 1
    elseif rone < (r[1]+r[2]+r[3]+r[4]+r[5]+r[6])
        a += 1
        B += 2
        updateB = 1
    elseif rone < (r[1]+r[2]+r[3]+r[4]+r[5]+r[6]+r[7])
        b -= 1
        A -= 2 # with cooperativity
        updateA = 1
    else
        b += 1
        A += 2
        updateA = 1
    end

    # update time
    rtwo = rand()  # second random number
    t = t - log(rtwo)/sum(rates)

    if updateA == 1
        Atemp[i] = A
        TtempA[i] = t
        i += 1
        Na += 1
    elseif updateB == 1
        Btemp[j] = B
        TtempB[j] = t
        j += 1
        Nb += 1
    end

    if i == batchsizeA + 1
        At = vcat(At, Atemp) # these two lines are probably the biggest use of time in this program
        TA = vcat(TA, TtempA)
        i = 1
        Atemp = fill(0, batchsizeA) # repreallocate
        TtempA = fill(0.0, batchsizeA)
    end

    if j == batchsizeB + 1
        Bt = vcat(Bt, Btemp) # these two lines are probably the biggest use of time in this program
        TB = vcat(TB, TtempB)
        j = 1
        Btemp = fill(0, batchsizeB) # repreallocate
        TtempB = fill(0.0, batchsizeB)
    end

end

# Need to now cat the final incomplete batch, will test before I do this though
Averytemp = fill(0, i-1)
Bverytemp = fill(0, j-1)
TAverytemp = fill(0.0, i-1)
TBverytemp = fill(0.0, j-1)
# Select the valued entries
for l = 1:(i-1)
    Averytemp[l] = Atemp[l]
    TAverytemp[l] = TtempA[l]
end
# Select the valued entries
for l = 1:(j-1)
    Bverytemp[l] = Btemp[l]
    TBverytemp[l] = TtempB[l]
end

# clear temps
clear!(:Atemp)
clear!(:Btemp)
clear!(:TtempA)
clear!(:TtempB)
gc()

# Then cat onto vector
At = vcat(At, Averytemp)
Bt = vcat(Bt, Bverytemp)
TA = vcat(TA, TAverytemp)
TB = vcat(TB, TBverytemp)

# clear temps
clear!(:Averytemp)
clear!(:Bverytemp)
clear!(:TAverytemp)
clear!(:TBverytemp)
gc()

## Now find probability distribution
maxA = convert(Int64, maximum(At)+1)
PA = fill(0.0, maxA)
maxB = convert(Int64, maximum(Bt)+1)
PB = fill(0.0, maxB)

for i = 2:length(At)
    PA[1+At[i-1]] = PA[1+At[i-1]] + (TA[i]-TA[i-1])/t;
end

for i = 2:length(Bt)
    PB[1+Bt[i-1]] = PB[1+Bt[i-1]] + (TB[i]-TB[i-1])/t;
end

aveA = 0
for i = 1:length(PA)
    aveA = aveA + (i-1)*PA[i]
end

aveB = 0
for i = 1:length(PB)
    aveB = aveB + (i-1)*PB[i]
end

print("Gillespie done!\n")

gr() # set plotting back end to plotlyjs()
# Plot
pone = plot(TA, At, xlabel = "Time", ylabel = "Solution", legend = false)
pone = plot!(pone, TB, Bt)
print("First plot done!\n")

# Clear up At & Bt
clear!(:At)
clear!(:Bt)
gc() # collect garbage

# 1st histograms
binsA = 0:(length(PA)-1)
ptwo = plot(binsA, Na*PA, xlabel = "A", ylabel = "Frequency", color = :blue, linetype=:bar, xlim = (-1,3*Ω), legend = false)
annotate!(1.5*Ω,0.5*maximum(Na*PA),text("<A>=$(aveA)",:left))
print("Second plot done!\n")

# 2nd histogram
binsB = 0:(length(PB)-1)
pthree = plot(binsB, Nb*PB, xlabel = "B", ylabel = "Frequency", color = :red, linetype=:bar, xlim = (-1,3*Ω), legend = false)
annotate!(1.5*Ω,0.5*maximum(Nb*PB),text("<B>=$(aveB)",:left))
print("Third plot done!\n")

plot(pone, ptwo, pthree, layout=(3,1))#pthree

print("Plots combined!\n")

savefig("../Results/SolutionvsTime.png")
print("Plot Saved!\n")

# print moments to screen
print(aveA)
print("\n")
print(aveB)
print("\n")

print("Finished!\n")
