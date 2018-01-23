#!/usr/bin/env julia
# Julia Script for Gillespie simulation of toy driven ecosystem, this now
# only calculates the histogram rather than trajectories to save computation time

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
f = 0.01 # Promoter switching
rr = 0.01

Na = 0 # counts of A and B, mainly to give an idea of statisical weight later on
Nb = 0

Ti = 0.0  # initial time
Tf = 50000 # end of simulation time in s

TtempA = [Ti; Ti]
TtempB = [Ti; Ti]
PA = [0.0]
PB = [0.0]

t = TtempA[1] # initialise
A = Ω
B = 0

a = 1
b = 1

# Main loop
while t <= Tf
    # set prior A and Bs
    Aprev = A
    Bprev = B

    # rates
    rates = [K*A, a*k, Q*B, b*q, f*B*(B-1)*a, (1-a)*rr, f*A*(A-1)*b, (1-b)*rr]
    r = rates/sum(rates)

    rone = rand() # first random number

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
        TtempA[2] = TtempA[1]
        TtempA[1] = t
        Na += 1
        if Aprev > length(PA) - 1
            PAnew = fill(0.0, Aprev+1)
            for i=1:length(PA)
                PAnew[i] = PA[i]
            end
            PA = PAnew
            PA[Aprev+1] += TtempA[1] - TtempA[2]
        else
            PA[Aprev+1] += TtempA[1] - TtempA[2]
        end
    elseif updateB == 1
        TtempB[2] = TtempB[1]
        TtempB[1] = t
        Nb += 1
        if Bprev > length(PB) - 1
            PBnew = fill(0.0, Bprev+1)
            for i=1:length(PB)
                PBnew[i] = PB[i]
            end
            PB = PBnew
            PB[Bprev+1] += TtempB[1] - TtempB[2]
        else
            PB[Bprev+1] += TtempB[1] - TtempB[2]
        end
    end
end

# Now rescale PB and PA by t
PA = PA/t
PB = PB/t

aveA = 0
for i = 1:length(PA)
    aveA = aveA + (i-1)*PA[i]
end

aveB = 0
for i = 1:length(PB)
    aveB = aveB + (i-1)*PB[i]
end

print("Gillespie done!\n")

gr() # set plotting back end to gr()
# 1st histograms
binsA = 0:(length(PA)-1)
pone = plot(binsA, Na*PA, xlabel = "A", ylabel = "Frequency", color = :blue, linetype=:bar, xlim = (-1,3*Ω), legend = false)
annotate!(1.5*Ω,0.5*maximum(Na*PA),text("<A>=$(aveA)",:left))
print("First plot done!\n")

# 2nd histogram
binsB = 0:(length(PB)-1)
ptwo = plot(binsB, Nb*PB, xlabel = "B", ylabel = "Frequency", color = :red, linetype=:bar, xlim = (-1,3*Ω), legend = false)
annotate!(1.5*Ω,0.5*maximum(Nb*PB),text("<B>=$(aveB)",:left))
print("Second plot done!\n")

plot(pone, ptwo, layout=(2,1))

print("Plots combined!\n")

savefig("../Results/SolutionvsTime.png")
print("Plot Saved!\n")

# print moments to screen
print(aveA)
print("\n")
print(aveB)
print("\n")

print("Finished!\n")
