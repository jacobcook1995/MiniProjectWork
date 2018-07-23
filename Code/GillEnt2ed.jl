#!/usr/bin/env julia
# GillEnt2.jl
# A script to perform a Gillespie simulation and then calculate the entropy production
# of the paths by direct computation of the reversed probabilities.
# The script will then calculate steady state entropy via the Shannon method for comparison
# This will first be done for the two species case

using Roots
using Plots
using StatsBase
import GR

# Now construct the three relevant vectors of equations
function f!(F::Array{Float64,1},x::Array{Float64,1},k::Float64,q::Float64,K::Float64,Q::Float64,r::Float64,f::Float64)
    F[1] = k*r/(r+f*x[2]*x[2]) - K*x[1]
    F[2] = q*r/(r+f*x[1]*x[1]) - Q*x[2]
    return F
end

# Then construct the necessary matrices
function D!(D::Array{Float64,2},x::Array{Float64,1},k::Float64,q::Float64,K::Float64,Q::Float64,r::Float64,f::Float64)
    D[1,1] = k*r/(r+f*x[2]*x[2]) + K*x[1]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = q*r/(r+f*x[1]*x[1]) + Q*x[2]
    return D
end

function p!(p::Array{Float64,1},x::Array{Float64,1},k::Float64,q::Float64,kmin::Float64,qmin::Float64,r::Float64,f::Float64)
    p[1] = k*r/(r + f*x[2]^2) - kmin*x[1]
    p[2] = q*r/(r + f*x[1]^2) - qmin*x[2]
    return(p)
end

function d!(d::Array{Float64,1},x::Array{Float64,1},K::Float64,Kmin::Float64,Q::Float64,Qmin::Float64)
    d[1] = K*x[1] - Kmin
    d[2] = Q*x[2] - Qmin
    return(d)
end

# finction to find start end and saddle points
function nullcline(r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64)
    a = 2
    b = 2
    A1(x) = k*r/(K*(r+f*x^a))
    A2(x) = (r/f*(q/(Q*x)-1))^(1/b)
    g(x) = k*r/(K*(r+f*x^a)) - (r/f*(q/(Q*x)-1))^(1/b) #A1(x) - A2(x)
    xs = fzeros(g, 0, q/Q)
    ss1 = [A1(xs[1]); xs[1]]
    sad = [A1(xs[2]); xs[2]]
    ss2 = [A1(xs[3]); xs[3]]
    print(ss1)
    print("\n")
    print(sad)
    print("\n")
    print(ss2)
    print("\n")
    return (ss1,sad,ss2)
end

# function to construct the rates
function rates(A::Int64,B::Int64,k::Float64,K::Float64,q::Float64,Q::Float64,kmin::Float64,
                Kmin::Float64,qmin::Float64,Qmin::Float64,r::Float64,f::Float64)
    rates = [ r*k/(r + f*B*(B-1)), kmin*A, K*A, Kmin, r*q/(r + f*A*(A-1)), qmin*B, Q*B, Qmin ]
    return(rates)
end

# function to calculate the time step
function timstep(rates::AbstractVector)
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::AbstractVector,vars::Array{Int64,1},reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars[1] += 1 # A produced
        p = rs[1]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars[1] -= 1 # A unravels
        p = rs[2]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars[1] -= 1 # A decays
        p = rs[3]
        reac = 3
    elseif r < rs[1] + rs[2] + rs[3] + rs[4]
        vars[1] += 1 # A regenerated
        p = rs[4]
        reac = 4
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5]
        vars[2] += 1 # B produced
        p = rs[5]
        reac = 5
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6]
        vars[2] -= 1 # B unravels
        p = rs[6]
        reac = 6
    elseif r < rs[1] + rs[2] + rs[3] + rs[4] + rs[5] + rs[6] + rs[7]
        vars[2] -= 1 # B decays
        p = rs[7]
        reac = 7
    else
        vars[2] += 1 # B regenerated
        p = rs[8]
        reac = 8
    end
    return(vars,p,reac)
end

# function to find reverse probability
function rev(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 9
        if reac % 2 == 0 # even rates backwards
            return(rs[reac-1])
        else # odd forward
            return(rs[reac+1])
        end
    else
        error("Invalid reaction code returned")
    end
end

function Gillespie!(stead::Array{Int64,1},K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,
                    qmin::Float64,f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},times::Array{Float64,1},vars::Array{Int64,2},
                    qs::Array{Float64,2},fs::Array{Float64,2},Ds::Array{Float64,3},ds::Array{Float64,2},ps::Array{Float64,2})
    # Preallocate for simulation
    times[1] = 0
    vars[:,1] = stead
    reac = 0
    for i = 1:noits
        # calculate rates
        rs = rates(vars[1,i],vars[2,i],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
        if i != 1
            # now use reac to calculate reverse rate
            pb[i-1] = rev(rs,reac)
        end
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[:,i+1], pf[i], reac = step(rs,vars[:,i],reac)
        qs[1,i] = (vars[1,i+1] - vars[1,i])/(τ)
        qs[2,i] = (vars[2,i+1] - vars[2,i])/(τ)
        posA = (vars[1,i+1] + vars[1,i])/(2)
        posB = (vars[2,i+1] + vars[2,i])/(2)
        fs[:,i] = f!(fs[:,i],[posA, posB],k,q,K,Q,r,f)
        ds[:,i] = d!(ds[:,i],[posA, posB],K,Kmin,Q,Qmin)
        ps[:,i] = p!(ps[:,i],[posA, posB],k,q,kmin,qmin,r,f)
        Ds[:,:,i] = D!(Ds[:,:,i],[posA, posB],k,q,K,Q,r,f)
        Ds[:,:,i] = Ds[:,:,i]/τ # applies here as we need a τ in each expression, D will be inverted later
        # final reverse rate
        if i == noits
            rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            pb[end] = rev(rs,reac)
        end
    end
    return(pf,pb,times,vars,qs,fs,Ds,ds,ps)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,r::Float64,f::Float64,K::Float64,Q::Float64,
                    k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64,
                    highA::Array{Int64,1},highB::Array{Int64,1})
    SA = zeros(noruns)
    SB = zeros(noruns)
    PA = ones(BigFloat,noruns)
    PB = ones(BigFloat,noruns)
    qqA = zeros(noruns)
    qfA = zeros(noruns)
    ffA = zeros(noruns)
    qqB = zeros(noruns)
    qfB = zeros(noruns)
    ffB = zeros(noruns)
    ddA = zeros(noruns)
    pdA = zeros(noruns)
    ppA = zeros(noruns)
    ddB = zeros(noruns)
    pdB = zeros(noruns)
    ppB = zeros(noruns)
    tA = zeros(noruns)
    tB = zeros(noruns)
    # preallocating arrays used inside function
    pfA = zeros(noits)
    pbA = zeros(noits)
    timesA = zeros(noits+1)
    varsA = fill(0,2,noits+1)
    pfB = zeros(noits)
    pbB = zeros(noits)
    timesB = zeros(noits+1)
    varsB = fill(0,2,noits+1)
    qA = zeros(2,noits)
    fA = zeros(2,noits)
    DA = zeros(2,2,noits)
    qB = zeros(2,noits)
    fB = zeros(2,noits)
    DB = zeros(2,2,noits)
    dA = zeros(2,noits)
    pA = zeros(2,noits)
    dB = zeros(2,noits)
    pB = zeros(2,noits)
    for j = 1:noruns
        pfA, pbA, timesA, varsA, qA, fA, DA, dA, pA = Gillespie!(highA,K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,pfA,pbA,timesA,varsA,qA,fA,DA,dA,pA)
        pfB, pbB, timesB, varsB, qB, fB, DB, dB, pB = Gillespie!(highB,K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,pfB,pbB,timesB,varsB,qB,fB,DB,dB,pB)
        # calculate total entropy production
        for i = 1:noits
            SA[j] += (log(pfA[i]) - log(pbA[i])) # do probability weighting at each step
            SB[j] += (log(pfB[i]) - log(pbB[i])) # does this work????????
            PA[j] *= pfA[i]
            PB[j] *= pfB[i]
            for l = 1:2
                qfA[j] += qA[l,i]*fA[l,i]/DA[l,l,i]
                qqA[j] += qA[l,i]*qA[l,i]/DA[l,l,i]
                ffA[j] += fA[l,i]*fA[l,i]/DA[l,l,i]
                ppA[j] += pA[l,i]*pA[l,i]/DA[l,l,i]
                pdA[j] += pA[l,i]*dA[l,i]/DA[l,l,i]
                ddA[j] += dA[l,i]*dA[l,i]/DA[l,l,i]
                qfB[j] += qB[l,i]*fB[l,i]/DB[l,l,i]
                qqB[j] += qB[l,i]*qB[l,i]/DB[l,l,i]
                ffB[j] += fB[l,i]*fB[l,i]/DB[l,l,i]
                ppB[j] += pB[l,i]*pB[l,i]/DB[l,l,i]
                pdB[j] += pB[l,i]*dB[l,i]/DB[l,l,i]
                ddB[j] += dB[l,i]*dB[l,i]/DB[l,l,i]
            end
        end

        # convert total entropy production to entropy production rate
        SA[j] = SA[j]/timesA[end]
        SB[j] = SB[j]/timesB[end]
        qfA[j] = qfA[j]/timesA[end]
        qqA[j] = qqA[j]/timesA[end]
        ffA[j] = ffA[j]/timesA[end]
        qfB[j] = qfB[j]/timesB[end]
        qqB[j] = qqB[j]/timesB[end]
        ffB[j] = ffB[j]/timesB[end]
        ppA[j] = ppA[j]/timesA[end]
        pdA[j] = pdA[j]/timesA[end]
        ddA[j] = ddA[j]/timesA[end]
        ppB[j] = ppB[j]/timesB[end]
        pdB[j] = pdB[j]/timesB[end]
        ddB[j] = ddB[j]/timesB[end]
        tA[j] = timesA[end]
        tB[j] = timesB[end]
    end
    println("Gillespies Done!")
    return(SA,SB,PA,PB,qfA,qqA,ffA,qfB,qqB,ffB,ppA,pdA,ddA,ppB,pdB,ddB,tA,tB)
end

# function to compute the shannon entropy at a fixed point
function shannon(point::Array{Float64,1},r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,
                kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    pa = k*r/(r+f*(point[2]-1)*point[2])
    pamin = kmin*point[1]
    pb = q*r/(r+f*(point[1]-1)*point[1])
    pbmin = qmin*point[2]
    da = K*point[1]
    damin = Kmin
    db = Q*point[2]
    dbmin = Qmin
    S1 = (pa - pamin)*log(pa/pamin)
    S2 = (pb - pbmin)*log(pb/pbmin)
    S3 = (da - damin)*log(da/damin)
    S4 = (db - dbmin)*log(db/dbmin)
    S = S1 + S2 + S3 + S4
    return(S)
end

# main function
function main()
    # General parameters
    Ω = 150.0
    K = 10.0
    k = K*Ω # steady state for A=k/K=1
    Q = 1.0
    q = Q*Ω
    kmin = 10.0^-20 # set all too 10.0^-20 for now
    Kmin = 10.0^-20
    qmin = 10.0^-20
    Qmin = 10.0^-20
    f = 1000.0/(Ω^2) # Promoter switching
    r = 10.0

    # first need to use these parameters to find a steady state
    star1, _, fin1 = nullcline(r,f,K,Q,k,q,kmin,qmin)
    # round star so that it becomes a vector of integers
    star2 = fill(0,2)
    fin2 = fill(0,2)
    for i = 1:2
        star2[i] = round(Int64,star1[i])
        fin2[i] = round(Int64,fin1[i])
    end
    # now the steady state will be used to find shannon entropy productions
    SAS = shannon(star1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    SBS = shannon(fin1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println("Entropy production rate of high A state via Shannon formula = $(SAS)")
    println("Entropy production rate of high B state via Shannon formula = $(SBS)")
    # now run multiple Gillespie simulations
    noits = 500000 # this number must be kept low to insure that the paths do not become less probable than the computer can ennumerate
    noruns = 2500#0 # memory allocation real problem if this is too large
    SA, SB, PA, PB, qfA, qqA, ffA, qfB, qqB, ffB, ppA, pdA, ddA, ppB, pdB, ddB, tA, tB = multgill(noits,noruns,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin,star2,fin2)
    # alter probabilies to match paths of equal length
    tAm = maximum(tA)
    tBm = maximum(tB)
    for i = 1:length(PA)
        PA[i] = (PA[i])^(tAm/tA[i])
    end
    for i = 1:length(PB)
        PB[i] = (PB[i])^(tBm/tB[i])
    end
    # renormalise probability
    PA = PA/sum(PA)
    PB = PB/sum(PB)
    SAG = SBG = qfAw = qfBw = ppAw = ppBw = ddAw = ddBw = pdAw = pdBw = qqAw = ffAw = qqBw = ffBw = 0
    for i = 1:noruns
        SAG += PA[i]*SA[i]
        SBG += PB[i]*SB[i]
        qfAw += PA[i]*qfA[i]
        qfBw += PB[i]*qfB[i]
        ppAw += PA[i]*ppA[i]
        ppBw += PB[i]*ppB[i]
        ddAw += PA[i]*ddA[i]
        ddBw += PB[i]*ddB[i]
        pdAw += PA[i]*pdA[i]
        pdBw += PB[i]*pdB[i]
        qqAw += PA[i]*qqA[i]
        ffAw += PA[i]*ffA[i]
        qqBw += PB[i]*qqB[i]
        ffBw += PB[i]*ffB[i]
    end
    # SA histogram
    pone = histogram(SA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,SAG)], color = :red)
    pone = vline!(pone, [sum(SA)/noruns], color = :green)
    savefig("../Results/HistSA$(ARGS[1]).png")

    # SB histogram
    pone = histogram(SB, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,SBG)], color = :red)
    pone = vline!(pone, [sum(SB)/noruns], color = :green)
    savefig("../Results/HistSB$(ARGS[1]).png")

    # qAf histogram
    pone = histogram(2*qfA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,2*qfAw)], color = :red)
    pone = vline!(pone, [2*sum(qfA)/noruns], color = :green)
    savefig("../Results/HistqfA$(ARGS[1]).png")

    # qBf histogram
    pone = histogram(2*qfB, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,2*qfBw)], color = :red)
    pone = vline!(pone, [2*sum(qfB)/noruns], color = :green)
    savefig("../Results/HistqfB$(ARGS[1]).png")

    # ppA ddA histogram
    epA = 2*(ppA+ddA)
    h = fit(Histogram, epA, closed=:right)
    pone = histogram(epA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,2*(ppAw+ddAw))], color = :red)
    pone = vline!(pone, [sum(epA)/noruns], color = :green)
    annotate!(minimum(epA),0.15*maximum(h.weights),text("mean ratio=$((sum(SA))/sum(epA))\nweighted mean ratio=$(convert(Float64,SAG)/convert(Float64,2*(ppAw+ddAw)))",:left))
    savefig("../Results/HistppAddA$(ARGS[1]).png")

    # ppA ddA histogram
    epB = 2*(ppB+ddB)
    h = fit(Histogram, epB, closed=:right)
    pone = histogram(epB, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,2*(ppBw+ddBw))], color = :red)
    pone = vline!(pone, [sum(epB)/noruns], color = :green)
    annotate!(minimum(epB),0.15*maximum(h.weights),text("mean ratio=$((sum(SB))/sum(epB))\nweighted mean ratio=$(convert(Float64,SBG)/convert(Float64,2*(ppBw+ddBw)))",:left))
    savefig("../Results/HistppBddB$(ARGS[1]).png")

    # pdA histogram
    pone = histogram(4*(pdA), xlabel = "Entropy Flow", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,4*(pdAw))], color = :red)
    pone = vline!(pone, [4*sum(pdA)/noruns], color = :green)
    savefig("../Results/HistpdA$(ARGS[1]).png")

    # pdB histogram
    pone = histogram(4*(pdB), xlabel = "Entropy Flow", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,4*(pdBw))], color = :red)
    pone = vline!(pone, [4*sum(pdB)/noruns], color = :green)
    savefig("../Results/HistpdB$(ARGS[1]).png")

    # action steady state A histogram
    ActA = 0.5*(qqA+ffA)-qfA
    h = fit(Histogram, ActA, closed=:right)
    pone = histogram(ActA, xlabel = "Action", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,0.5*(qqAw+ffAw) - qfAw)], color = :red)
    pone = vline!(pone, [(sum(ActA))/noruns], color = :green)
    println(maximum(h.weights))
    annotate!(0.4*maximum(ActA),0.05*maximum(h.weights),text("weighted mean Action=$(convert(Float64,0.5*(qqAw+ffAw) - qfAw))\nmean Action=$((sum(ActA))/noruns)",:left))
    savefig("../Results/HistActA$(ARGS[1])")

    # action steady state B histogram
    ActB = 0.5*(qqB+ffB)-qfB
    h = fit(Histogram, ActB, closed=:right)
    pone = histogram(ActB, xlabel = "Action", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [convert(Float64,0.5*(qqBw+ffBw) - qfBw)], color = :red)
    pone = vline!(pone, [(sum(ActB))/noruns], color = :green)
    println(maximum(h.weights))
    annotate!(0.4*maximum(ActB),0.05*maximum(h.weights),text("weighted mean Action=$(convert(Float64,0.5*(qqBw+ffBw) - qfBw))\nmean Action=$((sum(ActB))/noruns)",:left))
    savefig("../Results/HistActB$(ARGS[1])")
    return(nothing)
end

@time main()
