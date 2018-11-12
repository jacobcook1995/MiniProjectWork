#!/usr/bin/env julia
# GillSchlogl.jl
# A script to perform a Gillespie simulation for the Schlogl model

using Roots
using Plots
import GR

# Now construct the three relevant vectors of equations
function f!(F::Float64,X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    F = k1*V - K1*X - k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return F
end

# Now construct the three relevant vectors of equations
function f!(F::Float64,X::Float64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    F = k1*V - K1*X - k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return F
end

# Then construct the necessary matrices
function D!(D::Float64,X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    D = k1*V + K1*X + k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return D
end

# Then construct the necessary matrices
function D!(D::Float64,X::Float64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    D = k1*V + K1*X + k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return D
end

# finction to find start end and saddle points
function nullcline(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,high2low::Bool,V::Float64)
    # write out equation to be solved
    f(x) = k1 - K1*x - k2*(x^3) + K2*B*(x^2)
    three = false
    Xs = [ 0.0, 0.0, 0.0 ]
    while three == false
        X = fzeros(f, 0, 10, no_pts = 1000)
        if length(X) == 3
            three = true
            Xs = X
        end
    end
    Xs = Xs*V
    sad = Xs[2]
    if high2low == true
        ss1 = Xs[3]
        ss2 = Xs[1]
    else
        ss1 = Xs[1]
        ss2 = Xs[3]
    end
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

# function to construct the rates
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64)
    rates = [ k1*V, K1*X, k2*X*(X-1)*(X-2)/(V*V), K2*B*X*(X-1)/V]
    return(rates)
end

# function to calculate the time step
function timstep(rates::AbstractVector)
    r = rand()
    τ = -log(r)/sum(rates)
    return(τ)
end

# function to advance gillespie one step
function step(rates::AbstractVector,vars::Int64,reac::Int64)
    r = rand()
    rs = rates/sum(rates)
    p = 0 # probability used for forward path
    if r < rs[1]
        vars += 1 # X produced from A (not simulated)
        p = rs[1]
        reac = 1
    elseif r < rs[1] + rs[2]
        vars -= 1 # X unravels to A (not simulated)
        p = rs[2]
        reac = 2
    elseif r < rs[1] + rs[2] + rs[3]
        vars -= 1 # X decays to B (not simulated)
        p = rs[3]
        reac = 3
    else
        vars += 1 # X regenerated from B (not simulated)
        p = rs[4]
        reac = 4
    end
    return(vars,p,reac)
end

# function to find reverse probability
function rev(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 5
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
    if reac > 0 || reac < 5
        if reac == 2 || reac == 3
            return(rs[1]+rs[4])
        else
            return(rs[2]+rs[3])
        end
    else
        error("Invalid reaction code returned")
    end
end
# function to find forward probability from reduced master equation
function forr(rs::AbstractVector,reac::Int64)
    rs = rs/sum(rs)
    if reac > 0 || reac < 5
        if reac == 2 || reac == 3
            return(rs[2]+rs[3])
        else
            return(rs[1]+rs[4])
        end
    else
        error("Invalid reaction code returned")
    end
end

# function to calculate entropy production via Schnakenberg method for full master equation
function schf(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64,X::Float64)
    r1 = k1*V
    rmin1 = K1*X
    r2 = k2*(X^3)/(V^2)
    rmin2 = K2*B*(X^2)/V
    F1 = r1 - rmin1
    A1 = log(r1/rmin1)
    F2 = r2 - rmin2
    A2 = log(r2/rmin2)
    S1 = F1*A1
    S2 = F2*A2
    S = (S1 + S2)/V # volume rescale
    return(S)
end

# function to calculate entropy production via Schnakenberg method for reduced master equation
function schr(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Float64,X::Float64)
    r1 = k1*V + K2*B*(X^2)/V
    rmin1 = K1*X + k2*(X^3)/(V^2)
    F = r1 - rmin1
    A = log(r1/rmin1)
    S = F*A/V # volume rescale
    return(S)
end

function Gillespie!(stead::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},pfr::Array{Float64,1},pbr::Array{Float64,1},
                    times::Array{Float64,1},vars::Array{Int64,1},V::Float64,qs::Array{Float64,1},fs::Array{Float64,1},Ds::Array{Float64,1})
    # Preallocate for simulation
    times[1] = 0
    vars[1] = stead
    reac = 0
    tup = tdown = 0
    X = zeros(noits)
    for i = 1:noits
        # calculate rates
        rs = rates(vars[i],k1,K1,k2,K2,B,V)
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
        vars[i+1], pf[i], reac = step(rs,vars[i],reac)
        pfr[i] = forr(rs,reac)
        posX = (vars[i+1] + vars[i])/(2)#vars[i]
        qs[i] = (vars[i+1] - vars[i])/(V*τ)
        fs[i] = f!(fs[i],posX,k1,K1,k2,K2,B,V)
        Ds[i] = D!(Ds[i],posX,k1,K1,k2,K2,B,V)
        Ds[i] = Ds[i]/τ # applies here as we need a τ in each expression, D will be inverted later
        # final reverse rate
        if i == noits
            rs = rates(vars[end],k1,K1,k2,K2,B,V)
            pb[end] = rev(rs,reac)
            pbr[end] = revr(rs,reac)
        end
    end
    return(pf,pb,pfr,pbr,times,vars,qs,fs,Ds)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,
                    star1::Float64,mid1::Float64,fin1::Float64,V::Float64)
    SX = zeros(noruns)
    SXr = zeros(noruns)
    minX = zeros(noruns)
    maxX = zeros(noruns)
    qq = zeros(noruns)
    qf = zeros(noruns)
    ff = zeros(noruns)
    ts = zeros(noruns)
    p = zeros(BigFloat,noruns)
    qs = zeros(noits)
    fs = zeros(noits)
    Ds = zeros(noits)
    # generate high and low states
    star2 = round(Int64,star1)
    fin2 = round(Int64,fin1)
    mid2 = round(Int64,mid1)
    # preallocating arrays used inside function
    pfX = zeros(noits)
    pbX = zeros(noits)
    pfXr = zeros(noits)
    pbXr = zeros(noits)
    timesX = zeros(noits+1)
    varsX = fill(0,noits+1)
    for j = 1:noruns
        p[j] = 1
        init = star2 # generate an initial condition here, should be more complex in future
        pfX, pbX, pfXr, pbXr, timesX, varsX, qs, fs, Ds = Gillespie!(init,k1,K1,k2,K2,B,noits,pfX,pbX,pfXr,pbXr,timesX,varsX,V,qs,fs,Ds)
        # calculate total entropy production
        for i = 1:noits
            SX[j] += log(pfX[i]) - log(pbX[i])
            SXr[j] += log(pfXr[i]) - log(pbXr[i])
            p[j] *= pfX[i]
            qf[j] += qs[i]*fs[i]/Ds[i]
            qq[j] += qs[i]*qs[i]/Ds[i]
            ff[j] += fs[i]*fs[i]/Ds[i]
        end
        # convert total entropy production to entropy production rate
        SX[j] = SX[j]/timesX[end]
        SXr[j] = SXr[j]/timesX[end]
        ts[j] = timesX[end]
        qf[j] = qf[j]/timesX[end]
        qq[j] = qq[j]/timesX[end]
        ff[j] = ff[j]/timesX[end]
        # store this data so validity can be checked later on
        minX[j] = minimum(varsX)
        maxX[j] = maximum(varsX)
    end
    println("Gillespies Done!")
    return(SX,SXr,minX,maxX,p,ts,qf,qq,ff)
end

# main function
function main()
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    B = 4.0
    V = 2500.0
    high2low = false
    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(k1,K1,k2,K2,B,high2low,V)
    # now run multiple Gillespie simulations
    noits = 2500#00
    noruns = 500
    SX, SXr, minX, maxX, pf, ts, qf, qq, ff = multgill(noits,noruns,k1,K1,k2,K2,B,star1,mid1,fin1,V)
    scatter(SXr/V,2*qf)
    savefig("../Results/test.png")
    println(sum(SX)/(V*noruns))
    println(sum(SXr)/(V*noruns))
    println(2*sum(qf)/(noruns))
    # calculations here
    # rescale pf here
    t = maximum(ts)
    println(t)
    # qf histogram
    pone = histogram(2*qf, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [2*sum(qf)/noruns], color = :green)
    savefig("../Results/Histqf$(ARGS[1]).png")
    # ff histogram
    pone = histogram(ff, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [sum(ff)/noruns], color = :green)
    savefig("../Results/Histff$(ARGS[1]).png")
    # qq histogram
    pone = histogram(qq, xlabel = "Entropy Production", ylabel = "Frequency", color = :blue, legend = false)
    pone = vline!(pone, [sum(qq)/noruns], color = :green)
    savefig("../Results/Histqq$(ARGS[1]).png")
    return(nothing)
end

@time main()
