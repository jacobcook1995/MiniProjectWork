#!/usr/bin/env julia
# GillEnt2.jl
# A script to perform a Gillespie simulation and then calculate the entropy production
# of the paths by direct computation of the reversed probabilities.
# The script will then calculate steady state entropy via the Shannon method for comparison
# This will first be done for the two species case

using Roots
using Plots
using StatsBase
using SymPy
import GR

# Now construct the three relevant vectors of equations
function f!(F::Array{Float64,1},x::Array{Float64,1},k::Float64,q::Float64,K::Float64,Q::Float64,r::Float64,f::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    F[1] = k*r/(r+f*x[2]*x[2]) - (K+kmin)*x[1] + Kmin
    F[2] = q*r/(r+f*x[1]*x[1]) - (Q+qmin)*x[2] + Qmin
    return F
end

# Then construct the necessary matrices
function D!(D::Array{Float64,2},x::Array{Float64,1},k::Float64,q::Float64,K::Float64,Q::Float64,r::Float64,f::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    D[1,1] = k*r/(r+f*x[2]*x[2]) + (K+kmin)*x[1] + Kmin
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = q*r/(r+f*x[1]*x[1]) + (Q+qmin)*x[2] + Qmin
    return D
end

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function nullcline(ps::Array{Float64,1},high2low::Bool)
    # gonna use SymPy here
    # A1(x) = sqrt.((r/f)*(q./((qmin+Q).*x-Qmin) - 1))
    # A2(x) = (1/(kmin+K))*((k*r)./(r+f*x.^2) + Kmin)
    A1(x) = real(sqrt(complex((ps[10]/ps[9])*(ps[4]/((ps[3]+ps[7])*x - ps[8]) - 1))))
    A2(x) = (1/(ps[5]+ps[1]))*((ps[2]*ps[10])/(ps[10]+ps[9]*x^2) + ps[6])
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    while three == false
        bs1 = fzeros(g, 0.0, 0.1) # catches zeros about the origin
        bs2 = fzeros(g, 0.1, 2*ps[4]/ps[3]) # the upper bound here is potentially problematic
        bs = vcat(bs1,bs2)
        n = length(bs)
        gs = 0
        bad = zeros(Int64,0)
        for i = 1:n
            # check if this is actual solution and not artifact of the numerical method
            gs = g(bs[i])
            tol = 1.0e-14
            if gs >= 0 + tol || gs <= 0 - tol
                bad = append!(bad,i)
            end
        end
        if length(bad) != 0
            n = n - length(bad)
            bs = deleteat!(bs, bad)
        end
        if n == 3
            three = true
        end
    end
    sad = [ A1(bs[2]), bs[2] ]
    if high2low == true
        ss1 = [ A1(bs[1]), bs[1] ]
        ss2 = [ A1(bs[3]), bs[3] ]
    else
        ss1 = [ A1(bs[3]), bs[3] ]
        ss2 = [ A1(bs[1]), bs[1] ]
    end
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
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
# function to find reverse probability from reduced master equation
function revr(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 9
        if reac == 1 || reac == 4
            return(rs[2]+rs[3])
        elseif reac == 2 || reac == 3
            return(rs[1]+rs[4])
        elseif reac == 5 || reac == 8
            return(rs[6]+rs[7])
        else
            return(rs[5]+rs[8])
        end
    else
        error("Invalid reaction code returned")
    end
end
# function to find forward probability from reduced master equation
function forr(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 9
        if reac == 1 || reac == 4
            return(rs[1]+rs[4])
        elseif reac == 2 || reac == 3
            return(rs[2]+rs[3])
        elseif reac == 5 || reac == 8
            return(rs[5]+rs[8])
        else
            return(rs[6]+rs[7])
        end
    else
        error("Invalid reaction code returned")
    end
end

function Gillespie!(stead::Array{Int64,1},K::Float64,k::Float64,Q::Float64,q::Float64,kmin::Float64,
                    qmin::Float64,f::Float64,r::Float64,Kmin::Float64,Qmin::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},pfr::Array{Float64,1},pbr::Array{Float64,1},
                    times::Array{Float64,1},vars::Array{Int64,2},qs::Array{Float64,2},fs::Array{Float64,2},
                    Ds::Array{Float64,3},Ω::Float64)
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
            pbr[i-1] = revr(rs,reac)
        end
        # calculate timestep
        τ = timstep(rs)
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[:,i+1], pf[i], reac = step(rs,vars[:,i],reac)
        pfr[i] = forr(rs,reac)
        qs[1,i] = (vars[1,i+1] - vars[1,i])/(Ω*τ)
        qs[2,i] = (vars[2,i+1] - vars[2,i])/(Ω*τ)
        posA = (vars[1,i+1] + vars[1,i])/(2)
        posB = (vars[2,i+1] + vars[2,i])/(2)
        fs[:,i] = f!(fs[:,i],[posA,posB],k,q,K,Q,r,f,kmin,qmin,Kmin,Qmin)
        Ds[:,:,i] = D!(Ds[:,:,i],[posA,posB],k,q,K,Q,r,f,kmin,qmin,Kmin,Qmin)
        Ds[:,:,i] = Ds[:,:,i]/τ # applies here as we need a τ in each expression, D will be inverted later
        # final reverse rate
        if i == noits
            rs = rates(vars[1,end],vars[2,end],k,K,q,Q,kmin,Kmin,qmin,Qmin,r,f)
            pb[end] = rev(rs,reac)
            pbr[end] = revr(rs,reac)
        end
    end
    j = 0
    if vars[1,1] < vars[2,1]
        j = 1
    else
        j = 2
    end
    return(pf,pb,pfr,pbr,times,vars,qs,fs,Ds)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,r::Float64,f::Float64,K::Float64,Q::Float64,
                    k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64,
                    highA::Array{Int64,1},Ω::Float64)
    S = zeros(noruns)
    Sr = zeros(noruns)
    P = ones(BigFloat,noruns)
    qq = zeros(noruns)
    qf = zeros(noruns)
    ff = zeros(noruns)
    t = zeros(noruns)
    # preallocating arrays used inside function
    pf = zeros(noits)
    pbr = zeros(noits)
    pfr = zeros(noits)
    pb = zeros(noits)
    times = zeros(noits+1)
    vars = fill(0,2,noits+1)
    qs = zeros(2,noits)
    fs = zeros(2,noits)
    D = zeros(2,2,noits)

    for j = 1:noruns
        pf, pb, pfr, pbr, times, vars, qs, fs, D = Gillespie!(highA,K,k,Q,q,kmin,qmin,f,r,Kmin,Qmin,noits,pf,pb,pfr,pbr,times,vars,qs,fs,D,Ω)
        # calculate total entropy production
        for i = 1:noits
            S[j] += (log(pf[i]) - log(pb[i])) # do probability weighting at each step
            Sr[j] += (log(pfr[i]) - log(pbr[i]))
            P[j] *= pf[i]
            for l = 1:2
                qf[j] += qs[l,i]*fs[l,i]/D[l,l,i]
                qq[j] += qs[l,i]*qs[l,i]/D[l,l,i]
                ff[j] += fs[l,i]*fs[l,i]/D[l,l,i]
            end
        end

        # convert total entropy production to entropy production rate
        S[j] = S[j]/times[end]
        Sr[j] = Sr[j]/times[end]
        qf[j] = qf[j]/times[end]
        qq[j] = qq[j]/times[end]
        ff[j] = ff[j]/times[end]
        t[j] = times[end]
    end
    println("Gillespies Done!")
    return(S,Sr,P,qf,qq,ff,t)
end

# function to compute the schnakenberg entropy at a fixed point
function sch(point::Array{Float64,1},r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,q::Float64,
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
    Ω = 300.0
    K = 10.0
    k = K*Ω # steady state for A=k/K=1
    Q = 1.0
    q = Q*Ω
    kmin = 10.0^-10 # set all too 10.0^-20 for now
    Kmin = (10.0^-10)*(Ω)
    qmin = 10.0^-10
    Qmin = (10.0^-10)*(Ω)
    f = 1000.0/(Ω^2) # Promoter switching
    r = 10.0
    high2low = false

    # first need to use these parameters to find a steady state
    ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
    star1, _, fin1 = nullcline(ps,high2low)
    # round star so that it becomes a vector of integers
    star2 = fill(0,2)
    fin2 = fill(0,2)
    for i = 1:2
        star2[i] = round(Int64,star1[i])
        fin2[i] = round(Int64,fin1[i])
    end
    # now the steady state will be used to find schnakenberg entropy productions
    SAS = sch(star1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    SBS = sch(fin1,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println("Entropy production rate of high A state via Shannon formula = $(SAS)")
    println("Entropy production rate of high B state via Shannon formula = $(SBS)")
    # now run multiple Gillespie simulations
    noits = 100000#0 # this number must be kept low to insure that the paths do not become less probable than the computer can ennumerate
    noruns = 250#00 # memory allocation real problem if this is too large
    SA, SrA, PA, qfA, qqA, ffA, tA = multgill(noits,noruns,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin,star2,Ω)
    # alter probabilies to match paths of equal length
    tAm = maximum(tA)
    # renormalise probability
    PA = PA/sum(PA)
    # confirm matches schnakenberg
    println(sum(SA)/noruns)
    println(sum(SrA)/(noruns*Ω))
    println(2*sum(qfA)/noruns)

    # qAf histogram
    pone = histogram(2*qfA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [2*sum(qfA)/noruns], color = :green)
    savefig("../Results/HistqfA$(ARGS[1]).png")

    # reduced S histogram
    pone = histogram(SrA/Ω, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [sum(SrA)/(noruns*Ω)], color = :green)
    savefig("../Results/HistSrA$(ARGS[1]).png")

    # qqA histogram
    pone = histogram(qqA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [sum(qqA)/noruns], color = :green)
    savefig("../Results/HistqqA$(ARGS[1]).png")

    # ffA histogram
    pone = histogram(ffA, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [sum(ffA)/noruns], color = :green)
    savefig("../Results/HistffA$(ARGS[1]).png")

    return(nothing)
end

@time main()
