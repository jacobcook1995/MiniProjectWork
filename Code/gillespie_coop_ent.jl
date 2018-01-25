#!/usr/bin/env julia
# gillespie_coop_ent.jl
# Script that calculates the entropy production of the bistable gene switch

# Add package for plotting
using Plots
using Roots

# Parameters
const Ω = 30 # system size
const k = 100 # steady state for A=k/K=1
const K = k/Ω # K=k'
const q = 10 # steady state for B=q/Q=1
const Q = q/Ω # Q=q'
const f = 0.01 # Promoter switching
const r = 0.01

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
    entros = zeros(3,2) # gonna neglect switching at first
    decA = K*point[1]
    decArev = 10.0^-60 # Any protein is equally (very un-) likely to spontanously assemble from celluar enviroment
    decB = Q*point[2]
    decBrev = 10.0^-60
    prodA = k*r/(r+f*point[2]*point[2]) # compound rate based on adiabatic assumption
    prodArev = (10.0^-20)*point[1] # spontanous decay to unfolded state seems more likely
    prodB = q*r/(r+f*point[1]*point[1])
    prodBrev = (10.0^-20)*point[2]
    # Big question is what to do about a and b's here
    aon = (1 - a)*r + 10.0^-60
    aoff = f*point[2]*(point[2])*a + 10.0^-60# - 1)*a
    bon = (1 - b)*r + 10.0^-60
    boff = f*point[1]*(point[1])*b + 10.0^-60# - 1)*b
    entros[1,1] = (decA - decArev)*log(decA/decArev)
    entros[1,2] = (decB - decBrev)*log(decB/decBrev)
    entros[2,1] = (prodA - prodArev)*log(prodA/prodArev)
    entros[2,2] = (prodB - prodBrev)*log(prodB/prodBrev)
    # look at calculating the entropy production from promotor switching
    entros[3,1] = (aon - aoff)*log(aon/aoff)
    entros[3,2] = (bon - boff)*log(bon/boff)
    return(entros)
end

# run gillespie simulation from Ti to Tf
function gille(Ti,Tf)
    batchsizeA = 1000000
    batchsizeB = 1000000
    batchsizea = 1000000
    batchsizeb = 1000000

    Na = 0 # counts of A and B, mainly to give an idea of statisical weight later on
    Nb = 0
    na = 0
    nb = 0

    At = []
    Bt = []
    at = []
    bt = []
    TA = []
    TB = []
    Ta = []
    Tb = []

    Atemp = fill(0, batchsizeA) # preallocate
    Btemp = fill(0, batchsizeB)
    atemp = fill(0, batchsizea)
    btemp = fill(0, batchsizeb)
    TtempA = fill(0.0, batchsizeA)
    TtempB = fill(0.0, batchsizeB)
    Ttempa = fill(0.0, batchsizea)
    Ttempb = fill(0.0, batchsizeb)

    Atemp[1] = Ω
    Btemp[1] = 0
    atemp[1] = 1
    btemp[1] = 1
    TtempA[1] = Ti
    TtempB[1] = Ti
    Ttempa[1] = Ti
    Ttempb[1] = Ti

    t = TtempA[1] # initialise
    A = Atemp[1]
    B = Btemp[1]
    a = atemp[1]
    b = btemp[1]

    i = 2 # start counts at two as the first element to fill is two
    j = 2
    k = 2
    l = 2

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
            atemp[k] = a
            Ttempa[k] = t
            k += 1
            na += 1
        elseif updateb == 1
            btemp[l] = b
            Ttempb[l] = t
            l += 1
            nb += 1
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

        if k == batchsizea + 1
            at = vcat(at, atemp)
            TA = vcat(Ta, Ttempa)
            k = 1
            atemp = fill(0, batchsizea)
            Ttempa = fill(0.0, batchsizea)
        end

        if l == batchsizeb + 1
            bt = vcat(bt, btemp)
            TB = vcat(Tb, Ttempb)
            l = 1
            btemp = fill(0, batchsizeb)
            Ttempb = fill(0.0, batchsizeb)
        end

    end

    # Need to now cat the final incomplete batch
    Averytemp = fill(0, i-1)
    Bverytemp = fill(0, j-1)
    averytemp = fill(0, k-1)
    bverytemp = fill(0, l-1)
    TAverytemp = fill(0.0, i-1)
    TBverytemp = fill(0.0, j-1)
    Taverytemp = fill(0.0, k-1)
    Tbverytemp = fill(0.0, l-1)
    # Select the valued entries
    for m = 1:(i-1)
        Averytemp[m] = Atemp[m]
        TAverytemp[m] = TtempA[m]
    end
    # Select the valued entries
    for m = 1:(j-1)
        Bverytemp[m] = Btemp[m]
        TBverytemp[m] = TtempB[m]
    end
    # Select the valued entries
    for m = 1:(k-1)
        averytemp[m] = atemp[m]
        Taverytemp[m] = Ttempa[m]
    end
    # Select the valued entries
    for m = 1:(l-1)
        bverytemp[m] = btemp[m]
        Tbverytemp[m] = Ttempb[m]
    end

    # clear temps
    clear!(:Atemp)
    clear!(:Btemp)
    clear!(:atemp)
    clear!(:btemp)
    clear!(:TtempA)
    clear!(:TtempB)
    clear!(:Ttempa)
    clear!(:Ttempb)
    gc()

    # Then cat onto vector
    At = vcat(At, Averytemp)
    Bt = vcat(Bt, Bverytemp)
    at = vcat(at, averytemp)
    bt = vcat(bt, averytemp)
    TA = vcat(TA, TAverytemp)
    TB = vcat(TB, TBverytemp)
    Ta = vcat(Ta, Taverytemp)
    Tb = vcat(Tb, Tbverytemp)

    # clear temps
    clear!(:Averytemp)
    clear!(:Bverytemp)
    clear!(:averytemp)
    clear!(:bverytemp)
    clear!(:TAverytemp)
    clear!(:TBverytemp)
    clear!(:Taverytemp)
    clear!(:Tbverytemp)
    gc()

    ## Now find probability distribution
    maxA = convert(Int64, maximum(At)+1)
    PA = fill(0.0, maxA)
    maxB = convert(Int64, maximum(Bt)+1)
    PB = fill(0.0, maxB)
    Pa = fill(0.0, 2)
    Pb = fill(0.0, 2)

    for i = 2:length(At)
        PA[1+At[i-1]] = PA[1+At[i-1]] + (TA[i]-TA[i-1])/t;
    end

    for i = 2:length(Bt)
        PB[1+Bt[i-1]] = PB[1+Bt[i-1]] + (TB[i]-TB[i-1])/t;
    end

    for i = 2:length(at)
        Pa[1+at[i-1]] = Pa[1+at[i-1]] + (Ta[i]-Ta[i-1])/t;
    end

    for i = 2:length(bt)
        Pb[1+bt[i-1]] = Pb[1+bt[i-1]] + (Tb[i]-Tb[i-1])/t;
    end

    aveA = 0
    for i = 1:length(PA)
        aveA = aveA + (i-1)*PA[i]
    end

    aveB = 0
    for i = 1:length(PB)
        aveB = aveB + (i-1)*PB[i]
    end

    avea = 0
    for i = 1:length(Pa)
        avea = avea + (i-1)*Pa[i]
    end

    aveb = 0
    for i = 1:length(Pb)
        aveb = aveb + (i-1)*Pb[i]
    end

    print("Gillespie done!\n")
    return(aveA, aveB, avea, aveb)
end

@time A, B, a, b = gille(0.0, 100000)
print("$(A),$(B),$(a),$(b)\n")
# ss1, sad, ss2 = nullcline() # ss1 is the more driven steady state
