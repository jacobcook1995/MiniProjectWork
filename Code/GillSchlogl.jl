#!/usr/bin/env julia
# GillSchlogl.jl
# A script to perform a Gillespie simulation for the Schlogl model

using Roots
using Plots
import GR

# Now construct the three relevant vectors of equations
function f!(F::Float64,x::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    F = k1*V - K1*x - k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return F
end

# Then construct the necessary matrices
function D!(D::Float64,x::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    D = k1*V + K1*x + k2*X*(X-1)*(X-2)/(V*V) + K2*B*X*(X-1)/V
    return D
end

# finction to find start end and saddle points
function nullcline(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,high2low::Bool,V::Int64)
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
function rates(X::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
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

function Gillespie!(stead::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,noits::Int64,
                    pf::Array{Float64,1},pb::Array{Float64,1},times::Array{Float64,1},vars::Array{Int64,1},
                    V::Int64,up::Int64,down::Int64,qs::Array{Float64,1},fs::Array{Float64,1},Ds::Array{Float64,1})
    # Preallocate for simulation
    times[1] = 0
    vars[1] = stead
    reac = 0
    tup = tdown = 0
    hl = zeros(noits)
    for i = 1:noits
        # calculate rates
        rs = rates(vars[i],k1,K1,k2,K2,B,V)
        if i != 1
            # now use reac to calculate reverse rate
            pb[i-1] = rev(rs,reac)
        end
        # calculate timestep
        τ = timstep(rs)
        # update marginal times
        if vars[i] <= down
            tdown += τ
            hl[i] = 2
        elseif vars[i] >= up
            tup += τ
            hl[i] = 1
        end
        # update time
        times[i+1] = times[i] + τ
        # do gillepsie step
        vars[i+1], pf[i], reac = step(rs,vars[i],reac)
        posX = (vars[i+1] + vars[i])/(2)
        qs[i] = (vars[i+1] - vars[i])/(τ)
        fs[i] = f!(fs[i],posX,k1,K1,k2,K2,B,V)
        Ds[i] = D!(Ds[i],posX,k1,K1,k2,K2,B,V)
        Ds[i] = Ds[i]/τ # applies here as we need a τ in each expression, D will be inverted later
        # final reverse rate
        if i == noits
            rs = rates(vars[end],k1,K1,k2,K2,B,V)
            pb[end] = rev(rs,reac)
        end
    end
    pup = tup/times[end]
    pdown = tdown/times[end]
    return(pf,pb,times,vars,pup,pdown,hl,qs,fs,Ds)
end

# function to run multiple short gillespie simulations in order to improve sampling statistics
function multgill(noits::Int64,noruns::Int64,k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,
                    star1::Float64,mid1::Float64,fin1::Float64,V::Int64)
    SX = zeros(noruns)
    minX = zeros(noruns)
    maxX = zeros(noruns)
    pup = zeros(noruns)
    pdown = zeros(noruns)
    qq = zeros(noruns)
    qf = zeros(noruns)
    ff = zeros(noruns)
    Sup = zeros(noruns)
    Sdown = zeros(noruns)
    ts = zeros(noruns)
    p = zeros(BigFloat,noruns)
    qs = zeros(noits)
    fs = zeros(noits)
    Ds = zeros(noits)
    # generate high and low states
    star2 = round(Int64,star1)
    fin2 = round(Int64,fin1)
    mid2 = round(Int64,mid1)
    # upper and lower bounds for the threshold
    up = ceil(Int64,mid1)
    down = floor(Int64,mid1)
    up = down # set this for now DELETE later
    # preallocating arrays used inside function
    pfX = zeros(noits)
    pbX = zeros(noits)
    timesX = zeros(noits+1)
    varsX = fill(0,noits+1)
    for j = 1:noruns
        p[j] = 1
        init = mid2 # generate an initial condition here, should be more complex in future
        pfX, pbX, timesX, varsX, pup[j], pdown[j], hl, qs, fs, Ds = Gillespie!(init,k1,K1,k2,K2,B,noits,pfX,pbX,timesX,varsX,V,up,down,qs,fs,Ds)
        # calculate total entropy production
        for i = 1:noits
            SX[j] += log(pfX[i]) - log(pbX[i])
            p[j] *= pfX[i]
            if hl[i] == 1
                Sup[j] += log(pfX[i]) - log(pbX[i])
            elseif hl[i] == 2
                Sdown[j] += log(pfX[i]) - log(pbX[i])
            end
            qf[j] += qs[i]*fs[i]/Ds[i]
            qq[j] += qs[i]*qs[i]/Ds[i]
            ff[j] += fs[i]*fs[i]/Ds[i]
        end

        # convert total entropy production to entropy production rate
        SX[j] = SX[j]/timesX[end]
        Sup[j] = Sup[j]/timesX[end]
        Sdown[j] = Sdown[j]/timesX[end]
        ts[j] = timesX[end]
        qf[j] = qf[j]/timesX[end]
        qq[j] = qq[j]/timesX[end]
        ff[j] = ff[j]/timesX[end]
        # store this data so validity can be checked later on
        minX[j] = minimum(varsX)
        maxX[j] = maximum(varsX)
    end
    println("Gillespies Done!")
    return(SX,minX,maxX,pup,pdown,Sup,Sdown,p,ts,qf,qq,ff)
end

# main function
function main()
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    B = 4.0
    V = 15
    high2low = false
    # first need to use these parameters to find a steady state
    star1, mid1, fin1 = nullcline(k1,K1,k2,K2,B,high2low,V)
    # now run multiple Gillespie simulations
    noits = 250000
    noruns = 50
    SX, minX, maxX, pup, pdown, Sup, Sdown, pf, ts, qf, qq, ff = multgill(noits,noruns,k1,K1,k2,K2,B,star1,mid1,fin1,V)
    println(sum(SX)/(V*noruns))
    # calculations here
    # rescale pf here
    t = maximum(ts)
    println(t)
    for i = 1:noruns
        pf[i] = (pf[i])^(t/ts[i])
    end
    pone = scatter(exp.((Sup-Sdown)/V),pup./pdown)
    savefig("../Results/ScatterSch.png")
    pone = scatter(log.(pf),pup./pdown)
    savefig("../Results/Scatter2Sch.png")
    pone = scatter(log.(pf),exp.((Sup-Sdown)/V))
    savefig("../Results/Scatter3Sch.png")
    return(nothing)
end

@time main()
