#!/usr/bin/env julia
# Schlögl.jl
# script to calculate the entropy production from the Schlögl model by 4 different
# methods

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
function stocAct(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,Xs::Float64)
    Nmax = 10000
    T = zeros(Nmax+1)
    X = zeros(Nmax+1)
    F1 = zeros(Nmax)
    F2 = zeros(Nmax)
    D = zeros(Nmax)
    x = round(Xs*V)
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
        elseif Rn[1] <= Rand1 && Rand1 < Rn[1] + Rn[2]
            x -= 1
        elseif Rn[1] + Rn[2] <= Rand1 && Rand1 < Rn[1] + Rn[2] + Rn[3]
            x -= 1
        else
            x += 1
        end
        # add updated x to vector
        X[i+1] = x

        # time increase (2. random number)
        Rand2 = rand()
        dt = -log(Rand2)/Rtot;
        t += dt
        T[i+1] = t
    end
    Si = Se = 0
    # now calculate term in action
    for i = 1:Nmax
        Si += 2*(F1[i]^2 + F2[i]^2)*(T[i+1]-T[i])/D[i]
        Se += 4*(F1[i]*F2[i])*(T[i+1]-T[i])/D[i]
    end
    Si = Si/T[end]
    Se = Se/T[end]
    println("Max X = $(maximum(X))")
    println("Min X = $(minimum(X))")
    return(X,T,Si,Se)
end

# method 3 master equation
function master(k1A::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64)
    # max. number of partiles
    N = V*10
    P = ones(BigFloat,N) # probability

    # find probability distribution
    for i = 1:N
        p = 1
        for j = 1:i
            Wp = k1A*V + K2*B*(j)*(j-1)/V # W_+(j-1)
            Wm = K1*j + k2*(j)*(j-1)*(j-2)/V^2 # W_-(j)
            p = p*Wp/Wm # formula found by iteration of master eq at NESS
        end
        P[i] = p
    end
    # normalize
    P = P/sum(P)
    # maximum 1
    Pmax1 = 0
    imax1 = 1
    for i = 2:N
        if P[i] >= P[i-1]
            Pmax1 = P[i]
            imax1 = i
        else
            break
        end
    end
    # minimum
    Pmin = P[imax1]
    imin = imax1
    for i = imin+1:N
        if P[i] <= P[i-1]
           Pmin = P[i]
           imin = i
        else
            break
        end
    end
    # maximum 2
    Pmax2 = P[imin]
    imax2 = imin

    for i = imin+1:N
        if P[i] >= P[i-1]
            Pmax2 = P[i]
            imax2 = i
        else
            break
        end
    end
    return(P,imax1,imin,imax2)
end

function entp(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,P::Array{BigFloat,1})
    # entropy production
    s = 0
    N = length(P)
    for i = 1:N-1
        if i > 3
            Wp1 = k1*V
            Wm1 = K1*i
            Wp2 = k2*(i-1)*(i-2)*(i-3)/(V*V)
            Wm2 = K2*B*(i-1)*(i-2)/V

            s += 0.5*(Wp1*P[i-1]-Wm1*P[i])*log((Wp1*P[i-1])/(Wm1*P[i]))
            s += 0.5*(Wm2*P[i-1]-Wp2*P[i])*log((Wm2*P[i-1])/(Wp2*P[i]))

            Wp1 = k1*V
            Wm1 = K1*(i+1)
            Wp2 = k2*(i+1)*i*(i-1)/(V*V)
            Wm2 = K2*B*i*(i-1)/V

            s += 0.5*(Wm1*P[i+1]-Wp1*P[i])*log((Wm1*P[i+1])/(Wp1*P[i]))
            s += 0.5*(Wp2*P[i+1]-Wm2*P[i])*log((Wp2*P[i+1])/(Wm2*P[i]))
        elseif i == 3
            Wp1 = k1*V
            Wm1 = K1*i

            s += 0.5*(Wp1*P[i-1]-Wm1*P[i])*log((Wp1*P[i-1])/(Wm1*P[i]))

            Wp1 = k1*V
            Wm1 = K1*(i+1)
            Wp2 = k2*(i+1)*i*(i-1)/(V*V)
            Wm2 = K2*B*i*(i-1)/V

            s += 0.5*(Wm1*P[i+1]-Wp1*P[i])*log((Wm1*P[i+1])/(Wp1*P[i]))
            s += 0.5*(Wp2*P[i+1]-Wm2*P[i])*log((Wp2*P[i+1])/(Wm2*P[i]))

        elseif i == 2
            Wp1 = k1*V
            Wm1 = K1

            s += 0.5*(Wp1*P[i-1]-Wm1*P[i])*log((Wp1*P[i-1])/(Wm1*P[i]))

            Wp1 = k1*V
            Wm1 = K1*(i+1)

            s += 0.5*(Wm1*P[i+1]-Wp1*P[i])*log((Wm1*P[i+1])/(Wp1*P[i]))

        elseif i==1
            Wp1 = k1*V
            Wm1 = 0

            s += 0

            Wp1 = k1*V
            Wm1 = K1*(i+1)

            s += 0.5*(Wm1*P[i+1]-Wp1*P[i])*log((Wm1*P[i+1])/(Wp1*P[i]))
        end

    end
    s = convert(Float64,s)
    return(s)
end

# now the second method extracting entropy production for stochastic action
function ProbEnt(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,Xs::Float64)
    Nmax = 100000
    T = zeros(Nmax+1)
    X = zeros(Nmax+1)
    F1 = zeros(Nmax)
    F2 = zeros(Nmax)
    D = zeros(Nmax)
    Pf = ones(Nmax)
    Pb = ones(Nmax)
    x = round(Xs*V)
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
    S = 0
    for i = 1:Nmax
        S += log(Pf[i]) - log(Pb[i])
    end
    S = S/T[end]
    return(X,T,Pf,Pb,S)
end

function main()
    # parameters
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    # # run determenistic script
    X1, X2, X3 = mydet(k1,K1,k2,K2)
    println(X3)
    # # Now do stochastic action script
    B = 4.0
    V = 200
    X, T, Si, Se = stocAct(k1,K1,k2,K2,B,V,X3)
    println(Si/V)
    println(Se/V)
    println((Si-Se)/V)
    P, _, _, _ = master(k1,K1,k2,K2,B,V)
    S = entp(k1,K1,k2,K2,B,V,P)
    println(S/V)
    X, T, Pf, Pr, S = ProbEnt(k1,K1,k2,K2,B,V,X3)
    println(S/V)
    return(nothing)
end

@time main()
