#!/usr/bin/env julia
# gillespie_4baths.jl
# Julia Script for Gillespie simulation of toy driven ecosystem, this now
# only calculates the histogram rather than trajectories to save computation time
# In this case 4 resivours are considered, protiens A and B, waste W, and substrate S

# Add package for plotting
using Plots
import GR

# Parameters
const Ω = 30 # system size
const k = 100 # steady state for A=k/K=1
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const K = k/Ω # K=k'
const Kmin = 10.0^-20
const q = 10 # steady state for B=q/Q=1
const qmin = 10.0^-20
const Q = q/Ω # Q=q'
const Qmin = 10.0^-20
const f = 10 # Promoter switching
const r = 10
const F = 250 # removal rate

function gillespie()
    Na = 0 # counts of A and B, mainly to give an idea of statisical weight later on
    Nb = 0
    Nw = 0
    Ns = 0

    Ti = 0.0  # initial time
    Tf = 1000#000#0 # end of simulation time in s

    TtempA = [Ti; Ti]
    TtempB = [Ti; Ti]
    TtempW = [Ti; Ti]
    TtempS = [Ti; Ti]
    PA = [0.0]
    PB = [0.0]
    PW = [0.0]
    PS = [0.0]

    t = TtempA[1] # initialise
    A = 0
    B = 0
    W = 1000
    S = 1000

    # Main loop
    while t <= Tf
        # set prior A and Bs
        Aprev = A
        Bprev = B
        Wprev = W
        Sprev = S

        # rates
        rates = [r*k*S/(r+f*B*(B-1)), kmin*A, K*A, Kmin*W, r*q*S/(r+f*A*(A-1)), qmin*B, Q*B, Qmin*W, F]
        rs = rates/sum(rates)

        rone = rand() # first random number

        updateA = false
        updateB = false
        updateW = false
        updateS = false

        # which reaction?
        if rone < rs[1]
            updateA = true
            updateS = true
            A += 1
            S -= 1
        elseif rone < (rs[1]+rs[2])
            updateA = true
            updateS = true
            A -= 1
            S += 1
        elseif rone < (rs[1]+rs[2]+rs[3])
            updateA = true
            updateW = true
            A -= 1
            W += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4])
            updateA = true
            updateW = true
            A += 1
            W -= 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5])
            updateB = true
            updateS = true
            B += 1
            S -= 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6])
            updateB = true
            updateS = true
            B -= 1
            S += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7])
            updateB = true
            updateW = true
            B -= 1
            W += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8])
            updateB = true
            updateW = true
            B += 1
            W -= 1
        else
            updateS = true
            updateW = true
            S += 1
            W -= 1
        end

        # update time
        rtwo = rand()  # second random number
        t = t - log(rtwo)/sum(rates)

        if updateA == true
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
        elseif updateB == true
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

        # Same process for W and S
        if updateW == true
            TtempW[2] = TtempW[1]
            TtempW[1] = t
            Nw += 1
            if Wprev > length(PW) - 1
                PWnew = fill(0.0, Wprev+1)
                for i = 1:length(PW)
                    PWnew[i] = PW[i]
                end
                PW = PWnew
                PW[Wprev+1] += TtempW[1] - TtempW[2]
            else
                PW[Wprev+1] += TtempW[1] - TtempW[2]
            end
        elseif updateS == true
            TtempS[2] = TtempS[1]
            TtempS[1] = t
            Ns += 1
            if Sprev > length(PS) - 1
                PSnew = fill(0.0, Sprev+1)
                for i=1:length(PS)
                    PSnew[i] = PS[i]
                end
                PS = PSnew
                PS[Sprev+1] += TtempS[1] - TtempS[2]
            else
                PS[Sprev+1] += TtempS[1] - TtempS[2]
            end
        end
    end
    # Now rescale PB and PA by t
    PA = PA/t
    PB = PB/t
    PW = PW/t
    PS = PS/t

    aveA = 0
    for i = 1:length(PA)
        aveA = aveA + (i-1)*PA[i]
    end

    aveB = 0
    for i = 1:length(PB)
        aveB = aveB + (i-1)*PB[i]
    end

    aveW = 0
    for i = 1:length(PW)
        aveW = aveW + (i-1)*PW[i]
    end

    aveS = 0
    for i = 1:length(PS)
        aveS = aveS + (i-1)*PS[i]
    end

    print("Gillespie done!\n")

    # print moments to screen
    print("$(aveA)\n")
    print("$(aveB)\n")
    print("$(aveW)\n")
    print("$(aveS)\n")

    gr() # set plotting back end to gr()
    # 1st histograms
    binsA = 0:(length(PA)-1)
    pone = plot(binsA, Na*PA, xlabel = "A", ylabel = "Frequency", color = :blue, linetype=:bar, xlim = (-1,3*length(PA)), legend = false)
    annotate!(1.5*length(PA), 0.5*maximum(Na*PA), text("<A>=$(aveA)",:left))
    print("First plot done!\n")

    # 2nd histogram
    binsB = 0:(length(PB)-1)
    ptwo = plot(binsB, Nb*PB, xlabel = "B", ylabel = "Frequency", color = :red, linetype=:bar, xlim = (-1,3*length(PB)), legend = false)
    annotate!(1.5*length(PB), 0.5*maximum(Nb*PB), text("<B>=$(aveB)",:left))
    print("Second plot done!\n")

    binsW = 0:(length(PW)-1)
    pthree = plot(binsW, Nw*PW, xlabel = "W", ylabel = "Frequency", color = :brown, linetype=:bar, xlim = (-1,3*length(PW)), legend = false)
    annotate!(1.5*(length(PW)), 0.5*maximum(Nw*PW), text("<W>=$(aveW)",:left))
    print("Third plot done!\n")

    # 2nd histogram
    binsS = 0:(length(PS)-1)
    pfour = plot(binsS, Ns*PS, xlabel = "S", ylabel = "Frequency", color = :green, linetype=:bar, xlim = (-1,3*length(PS)), legend = false)
    annotate!(1.5*(length(PS)), 0.5*maximum(Ns*PS), text("<S>=$(aveS)",:left))
    print("Fourth plot done!\n")

    plot(pone, ptwo, pthree, pfour, layout=(4,1))

    print("Plots combined!\n")

    savefig("../Results/SolutionvsTime.png")
    print("Plot Saved!\n")

    print("Finished!\n")

end

# run me big function
@time gillespie()
