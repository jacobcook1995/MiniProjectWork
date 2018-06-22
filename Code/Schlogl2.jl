#!/usr/bin/env julia
# Schlögl2.jl
# script to calculate the entropy production from the Schlögl, creates scatter graphs
# of path probility vs entropy production for ft and action methods

using Roots
using Plots
import GR

# method 1 determenistic
function mydet(k1::Float64,K1::Float64,k2::Float64,K2::Float64)
    # make range of B to plot over
    Bmin = 0.0
    Bmax = 6.0
    dB = 0.05
    B = Bmin:dB:Bmax
    X1 = X2 = X3 = 0
    # write out equation to be solved
    f(x,Bs) = k1 - K1*x - k2*(x^3) + K2*Bs*(x^2)
    #pone = plot(bg=:white)
    for i = 1:length(B)
        # make a g(x) to be solved for each entry in B
        g(x) = f(x,B[i])
        X = fzeros(g, 0, 10, no_pts = 1000)
        for j = 1:length(X)
            if imag(X[j]) == 0
                if j == 1
                    col = :red
                elseif j == 2
                    col = :yellow
                else
                    col = :green
                end
                # entropy production
                Sdot = (k1 - K1*X[j])*log(k1/(K1*X[j]))
                Sdot += (k2*(X[j]^3) - K2*B[i]*(X[j]^2))*log((k2*X[j])/(K2*B[i]))
                if B[i] == 4.0
                    println(Sdot)
                    X1 = X[1]
                    X2 = X[2]
                    X3 = X[3]
                end
                #pone = scatter!(pone,[B[i]],[Sdot],yscale=:log10,fillto=1e-2,legend=false,color=col)
             end
        end
    end
    #savefig("../Results/SchloglDetEnt.png")
    return(X1,X2,X3)
end

# now the second method extracting entropy production for stochastic action
function ProbEnt(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,Xs::Float64)
    Nmax = 50000
    T = zeros(Nmax+1)
    X = zeros(Nmax+1)
    D = zeros(Nmax)
    Pf = ones(BigFloat,Nmax)
    Pb = ones(BigFloat,Nmax)
    F1 = zeros(Nmax)
    F2 = zeros(Nmax)
    D = zeros(Nmax)
    x = round(Xs*V)
    A = zeros(Nmax)
    t = 0
    X[1] = x
    T[1] = t
    for i = 1:Nmax
        # transition rates
        R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
        F1[i] = R[1] - R[2]
        F2[i] = R[3] - R[4]
        Rtot = sum(R)
        D[i] = Rtot
        Rn = R/Rtot
        # which one?
        Rand1 = rand()

        if 0 <= Rand1 && Rand1 < Rn[1]
            x += 1
            Pf[i] = Rn[1]
            R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
            Rn = R/sum(R)
            Pb[i] = Rn[2]
        elseif Rn[1] <= Rand1 && Rand1 < Rn[1] + Rn[2]
            x -= 1
            Pf[i] = Rn[2]
            R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
            Rn = R/sum(R)
            Pb[i] = Rn[1]
        elseif Rn[1]+Rn[2] <= Rand1 && Rand1 < Rn[1]+Rn[2]+Rn[3];
            x -= 1
            Pf[i] = Rn[3]
            R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
            Rn = R/sum(R)
            Pb[i] = Rn[4]
        else
            x += 1
            Pf[i] = Rn[4]
            R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
            Rn = R/sum(R)
            Pb[i] = Rn[3]
        end
        # add updated x to vector
        X[i+1] = x

        # time increase (2. random number)
        Rand2 = rand()
        dt = -log(Rand2)/Rtot;
        t += dt
        T[i+1] = t
    end
    # now use forward and backward probabilities to calculate entropy production
    S1 = S2 = S3 = 0
    for i = 1:Nmax
        S1 += log(Pf[i]) - log(Pb[i])
    end
    S1 = convert(Float64,S1)
    for i = 1:Nmax
        S2 += 2*(F1[i]^2 + F2[i]^2)*(T[i+1]-T[i])/D[i]
        S3 -= 4*(F1[i]*F2[i])*(T[i+1] - T[i])/D[i]
    end
    # S1 = S1/T[end]
    # S2 = S2/T[end]
    Pft = 1
    for i = 1:Nmax
        Pft *= Pf[i]
    end
    At = 0
    for i = 1:Nmax
        v = (X[i+1] - X[i])/(T[i+1] - T[i])
        At += (v - F1[i] - F2[i])*(v - F1[i] - F2[i])*(T[i+1] - T[i])/D[i]
    end
    #At = At/T[end]
    # not sure the above line is needed in terms of the action for state weighting
    return(X,Pft,S1,S2,S3,At)
end

# function to run multiple Gillespie simulations at a steady state then calculate
# path probabilities and entropy productions at each state
function scatters(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,Xs::Float64)
    N = 20#00 # number of points
    P = zeros(BigFloat,N)
    S1 = zeros(N)
    S2 = zeros(N)
    S3 = zeros(N)
    A = zeros(N)
    lnP = zeros(N)
    maxX = zeros(N)
    minX = zeros(N)
    for i = 1:N
        X, P[i], S1[i], S2[i], S3[i], A[i] = ProbEnt(k1,K1,k2,K2,B,V,Xs)
        lnP[i] = convert(Float64,log(P[i]))
        maxX[i] = maximum(X)
        minX[i] = minimum(X)
    end
    # scatter(P,S1,xaxis=:log10)
    # savefig("../Results/Scatter1.png")
    # scatter(P,S2,xaxis=:log10)
    # savefig("../Results/Scatter2.png")
    # scatter(P,S1,xaxis=:log10)
    # scatter!(P,S2,xaxis=:log10)
    # savefig("../Results/Scatter3.png")
    scatter(A,-lnP,xaxis=:log10)
    savefig("../Results/Scatter1p.png")
    scatter(A-S3,-lnP,xaxis=:log10)
    savefig("../Results/Scatter2p.png")
    scatter(A-S2,-lnP,xaxis=:log10)
    savefig("../Results/Scatter3p.png")
    scatter(S1,-lnP,xaxis=:log10)
    savefig("../Results/Scatter4p.png")
    # Entropy production vs probability
    scatter(S2,-lnP,xaxis=:log10)
    savefig("../Results/Scatter5p.png")
    # Entropy flow vs probability
    scatter(-S3,-lnP,xaxis=:log10)
    savefig("../Results/Scatter6p.png")
end

function main()
    # parameters
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    X1, X2, X3 = mydet(k1,K1,k2,K2)
    B = 4.0
    V = 200
    # now run multiple Gillespies
    scatters(k1,K1,k2,K2,B,V,X1)
    return(nothing)
end

@time main()
