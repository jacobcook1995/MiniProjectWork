#!/usr/bin/env julia
# gillespie_coop_ent.jl
# Script that calculates the entropy production of the bistable gene switch

# Add package for plotting
using Plots
using Roots
using Atom
import GR

# Parameters
const Ω = 30 # system size
const k = 100 # steady state for A=k/K=1
const K = k/Ω # K=k'
const q = 10 # steady state for B=q/Q=1
const Q = q/Ω # Q=q'
const f = 0.0001 # Promoter switching
const r = 0.0001

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
function nullcline()
    a = 2 # redefining a global constant here, but that's fine really
    b = 2
    A1(x) = k*r/(K*(r+f*x^a))
    A2(x) = (r/f*(q/(Q*x)-1))^(1/b)
    g(x) = k*r/(K*(r+f*x^a)) - (r/f*(q/(Q*x)-1))^(1/b) #A1(x) - A2(x)
    xs = fzeros(g, 0, q/Q)
    ss1 = [A1(xs[1]); xs[1]]
    sad = [A1(xs[2]); xs[2]]
    ss2 = [A1(xs[3]); xs[3]]
    return (ss1,sad,ss2)
end

# function to calculate the entropy production of a point in protein space
function entprods(point,a,b)
    entros = zeros(4) # gonna neglect switching at first
    decA = K*point[1] + 10.0^-60
    decArev = 10.0^-60 # Any protein is equally (very un-) likely to spontanously assemble from celluar enviroment
    decB = Q*point[2] + 10.0^-60
    decBrev = 10.0^-60
    prodA = k*a + 10.0^-60 # compound rate based on adiabatic assumption
    prodArev = (10.0^-30)*point[1] + 10.0^-60 # spontanous decay to unfolded state seems more likely
    prodB = q*b + 10.0^-60
    prodBrev = (10.0^-30)*point[2] + 10.0^-60
    entros[1] = (decA - decArev)*log(decA/decArev)
    entros[2] = (decB - decBrev)*log(decB/decBrev)
    entros[3] = (prodA - prodArev)*log(prodA/prodArev)
    entros[4] = (prodB - prodBrev)*log(prodB/prodBrev)
    # look at calculating the entropy production from promotor switching
    return(entros)
end

# run gillespie simulation from Ti to Tf
function gille(Ti,Tf)
    batchsizeA = 1000000
    batchsizeB = 1000000
    batchsizea = 1000000
    batchsizeb = 1000000
    batchsizeS = 4000000

    Na = 0 # counts of A and B, mainly to give an idea of statisical weight later on
    Nb = 0
    na = 0
    nb = 0

    At = []
    Bt = []
    at = []
    bt = []
    St = []
    TA = []
    TB = []
    Ta = []
    Tb = []
    TS = []

    Atemp = fill(0, batchsizeA) # preallocate
    Btemp = fill(0, batchsizeB)
    atemp = fill(0, batchsizea)
    btemp = fill(0, batchsizeb)
    Stemp = fill(0.0, (batchsizeS, 4))
    TtempA = fill(0.0, batchsizeA)
    TtempB = fill(0.0, batchsizeB)
    Ttempa = fill(0.0, batchsizea)
    Ttempb = fill(0.0, batchsizeb)
    TtempS = fill(0.0, batchsizeS)

    Atemp[1] = Ω
    Btemp[1] = 0
    atemp[1] = 1
    btemp[1] = 1
    Stemp[1,:] = entprods([Atemp[1] Btemp[1]], atemp[1], btemp[1]) # calculate initial values for entropy production
    TtempA[1] = Ti
    TtempB[1] = Ti
    Ttempa[1] = Ti
    Ttempb[1] = Ti
    TtempS[1] = Ti

    t = TtempA[1] # initialise
    A = Atemp[1]
    B = Btemp[1]
    a = atemp[1]
    b = btemp[1]

    i = 2 # start counts at two as the first element to fill is two
    j = 2
    l = 2
    m = 2
    n = 2

    # Main loop
    while t <= Tf
        # rates
        rates = [K*A, a*k, Q*B, b*q, f*B*(B-1)*a, (1-a)*r, f*A*(A-1)*b, (1-b)*r]
        rs = rates/sum(rates)

        rone = rand() #first random number

        updateA = false
        updateB = false
        updatea = false
        updateb = false

        # which reaction?
        if rone < rs[1]
            A -= 1
            updateA = true
        elseif rone < (rs[1]+rs[2])
            A += 1
            updateA = true
        elseif rone < (rs[1]+rs[2]+rs[3])
            B -= 1
            updateB = true
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4])
            B += 1
            updateB = true
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5])
            a -= 1
            B -= 2 # With cooperativity
            updateB = true
            updatea = true
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6])
            a += 1
            B += 2
            updateB = true
            updatea = true
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7])
            b -= 1
            A -= 2 # with cooperativity
            updateA = true
            updateb = true
        else
            b += 1
            A += 2
            updateA = true
            updateb = true
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

        if updatea == 1
            atemp[l] = a
            Ttempa[l] = t
            l += 1
            na += 1
        elseif updateb == 1
            btemp[m] = b
            Ttempb[m] = t
            m += 1
            nb += 1
        end

        # Entropy prod calculation step
        Stemp[n,:] = entprods([A B], a, b)
        TtempS[n] = t
        n += 1

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

        if l == batchsizea + 1
            at = vcat(at, atemp)
            Ta = vcat(Ta, Ttempa)
            l = 1
            atemp = fill(0, batchsizea)
            Ttempa = fill(0.0, batchsizea)
        end

        if m == batchsizeb + 1
            bt = vcat(bt, btemp)
            Tb = vcat(Tb, Ttempb)
            m = 1
            btemp = fill(0, batchsizeb)
            Ttempb = fill(0.0, batchsizeb)
        end

        if n == batchsizeS + 1
            St = vcat(St, Stemp)
            TS = vcat(TS, TtempS)
            n = 1
            Stemp = fill(0.0, (batchsizeS, 4))
            TtempS = fill(0.0, batchsizeS)
        end
    end

    # Need to now cat the final incomplete batch
    Averytemp = fill(0, i-1)
    Bverytemp = fill(0, j-1)
    averytemp = fill(0, l-1)
    bverytemp = fill(0, m-1)
    Sverytemp = fill(0.0, (n-1, 4))
    TAverytemp = fill(0.0, i-1)
    TBverytemp = fill(0.0, j-1)
    Taverytemp = fill(0.0, l-1)
    Tbverytemp = fill(0.0, m-1)
    TSverytemp = fill(0.0, n-1)
    # Select the valued entries
    for p = 1:(i-1)
        Averytemp[p] = Atemp[p]
        TAverytemp[p] = TtempA[p]
    end
    # Select the valued entries
    for p = 1:(j-1)
        Bverytemp[p] = Btemp[p]
        TBverytemp[p] = TtempB[p]
    end
    # Select the valued entries
    for p = 1:(l-1)
        averytemp[p] = atemp[p]
        Taverytemp[p] = Ttempa[p]
    end
    # Select the valued entries
    for p = 1:(m-1)
        bverytemp[p] = btemp[p]
        Tbverytemp[p] = Ttempb[p]
    end
    # Select the valued entries
    for p = 1:(n-1)
        Sverytemp[p,:] = Stemp[p,:]
        TSverytemp[p] = TtempS[p]
    end

    # clear temps
    clear!(:Atemp)
    clear!(:Btemp)
    clear!(:atemp)
    clear!(:btemp)
    clear!(:Stemp)
    clear!(:TtempA)
    clear!(:TtempB)
    clear!(:Ttempa)
    clear!(:Ttempb)
    clear!(:TtempS)
    gc()


    # Then cat onto vector
    At = vcat(At, Averytemp)
    Bt = vcat(Bt, Bverytemp)
    at = vcat(at, averytemp)
    bt = vcat(bt, bverytemp)
    St = vcat(St, Sverytemp)
    TA = vcat(TA, TAverytemp)
    TB = vcat(TB, TBverytemp)
    Ta = vcat(Ta, Taverytemp)
    Tb = vcat(Tb, Tbverytemp)
    TS = vcat(TS, TSverytemp)

    # clear temps
    clear!(:Averytemp)
    clear!(:Bverytemp)
    clear!(:averytemp)
    clear!(:bverytemp)
    clear!(:Sverytemp)
    clear!(:TAverytemp)
    clear!(:TBverytemp)
    clear!(:Taverytemp)
    clear!(:Tbverytemp)
    clear!(:TSverytemp)
    gc()

    ## Now find probability distribution
    maxA = convert(Int64, maximum(At)+1)
    PA = fill(0.0, maxA)
    maxB = convert(Int64, maximum(Bt)+1)
    PB = fill(0.0, maxB)
    Pa = fill(0.0, 2)
    Pb = fill(0.0, 2)

    # Calculate for each species
    for i = 2:length(At)
        PA[1+At[i-1]] += (TA[i]-TA[i-1])/t;
    end
    for i = 2:length(Bt)
        PB[1+Bt[i-1]] += (TB[i]-TB[i-1])/t;
    end
    for i = 2:length(at)
        Pa[1+at[i-1]] += (Ta[i]-Ta[i-1])/t;
    end

    for i = 2:length(bt)
        Pb[1+bt[i-1]] += (Tb[i]-Tb[i-1])/t;
    end
    aveA = 0
    for i = 1:length(PA)
        aveA += (i-1)*PA[i]
    end
    aveB = 0
    for i = 1:length(PB)
        aveB += (i-1)*PB[i]
    end
    avea = 0
    for i = 1:length(Pa)
        avea += (i-1)*Pa[i]
    end
    aveb = 0
    for i = 1:length(Pb)
        aveb += (i-1)*Pb[i]
    end

    print("Gillespie done!\n")

    averages = zeros(4)
    averages[1] = aveA
    averages[2] = aveB
    averages[3] = avea
    averages[4] = aveb

    gr()
    pone = plot(TS, St)
    savefig("../Results/EntropiesvsTime") # takes a really long time
    ptwo = plot(TA, At)
    pthree = plot(TB, Bt)
    pfour = plot(Ta, at)
    pfive = plot(Tb, bt)

    plot(ptwo, pthree, pfour, pfive, layout=(4,1))
    print("Plots combined!\n")
    savefig("../Results/SolutionvsTime.png") # savefig takes a really long time

    return(averages)
end

@time Averages  = gille(0.0, 1000000)
print("$(Averages[1]),$(Averages[2]),$(Averages[3]),$(Averages[4])\n")
# ss1, sad, ss2 = nullcline() # ss1 is the more driven steady state
