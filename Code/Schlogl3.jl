#!/usr/bin/env julia
# Schlögl3.jl
# script to calculate the entropy production from the Schlögl, now trying to calculate the
# higher order terms in the Scholgl model, so that I can compute the full action

using Roots
using Plots
using SymPy
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

function Fs()
    # k = k+1, kmin = k-1, K = k+2, Kmin = k-2
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    # Make a symbolic version of the matrix, needs no input in this case
    F = K*x^3 - Kmin*b*x^2 + kmin*x - k*a
    return(F)
end

function dFs()
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    F = Fs()
    dF = diff(F,x)
    return(dF)
end

function es()
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    e = sqrt(K*x^3 + Kmin*b*x^2 + kmin*x + k*a)
    return(e)
end

function des()
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    e = es()
    de = diff(e,x)
    return(de)
end

function inves()
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    e = es()
    inve = 1/e
    return(inve)
end

function Ds()
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    e = es()
    D = e*e
    return(D)
end

# ps are parameters
function equats(ps::Array{Float64,1})
    k, kmin, K, Kmin, a, b, x = symbols("k,kmin,K,Kmin,a,b,x")
    F = Fs()
    dF = dFs()
    e = es()
    de = des()
    inve = inves()
    D = Ds()
    F = subs(F,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    dF = subs(dF,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    e = subs(e,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    de = subs(de,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    inve = subs(inve,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    D = subs(D,k=>ps[1],kmin=>ps[2],K=>ps[3],Kmin=>ps[4],a=>ps[5],b=>ps[6])
    return(F,dF,e,de,inve,D)
end

function Act(F::Sym,dF::Sym,de::Sym,inve::Sym,D::Sym,X::Float64,V::Int64,q::Float64)
    x = symbols("x")
    F = subs(F,x=>X) |> float
    dF = subs(dF,x=>X) |> float
    de = subs(de,x=>X) |> float
    inve = subs(inve,x=>X) |> float
    D = subs(D,x=>X) |> float
    A1 = (V/2)*(q)*(q)/D
    A2 = (V/2)*F*F/D
    A3 = V*q*F/D
    return(A1,A2,A3)
end

# now the second method extracting entropy production for stochastic action
function stocAct(k1::Float64,K1::Float64,k2::Float64,K2::Float64,B::Float64,V::Int64,Xs::Float64,F::Sym,dF::Sym,de::Sym,inve::Sym,D::Sym)
    Nmax = 100000
    T = zeros(Nmax+1)
    X = zeros(Nmax+1)
    F1 = zeros(Nmax)
    F2 = zeros(Nmax)
    x = round(Xs*V) + 20
    A1 = zeros(Nmax)
    A2 = zeros(Nmax)
    A3 = zeros(Nmax)
    t = 0
    X[1] = x
    T[1] = t
    for i = 1:Nmax
        # transition rates
        R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
        F1[i] = R[1] - R[2]
        F2[i] = R[3] - R[4]
        Rtot = sum(R)
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
    return(X,T)
end

function main()
    k = 0.5 # k_{+1}[A]
    kmin = 3.0 # k_{-1}
    K = 1.0 # k_{+2}
    Kmin = 1.0 # k_{-2}
    b = 4.0
    a = 1.0
    V = 20000#0000
    # generate functions
    F, dF, e, de, inve, D = equats([k,kmin,K,Kmin,a,b])
    # find steady states
    X1, X2, X3 = mydet(k,kmin,K,Kmin)
    # both of these are defined in terms of concentrations
    println(X1)
    println(X2)
    println(X3)
    # Now run gillespie
    X, T = stocAct(k,kmin,K,Kmin,b,V,X3,F,dF,de,inve,D)
    # Now do analysis with the data binned
    # should now define a set time interval and use this too obtain the number of bins (round down)
    dT = 0.0005
    nobins = floor(Int64,T[end]/dT)
    println(nobins)
    A1 = zeros(nobins)
    A2 = zeros(nobins)
    A3 = zeros(nobins)
    poss = zeros(nobins)
    t = 0
    one = X[1] # set posistion to start from
    plot(X/V)
    savefig("../Results/Graph1.png")
    j = 1
    for i = 1:nobins
        t += dT
        # now need to find end posistion for bin
        found = false
        while found == false
            if T[j] < t && j != length(X)
                j += 1
            else
                found = true
            end
        end
        # linear interpolation step
        distf = T[j] - t
        distb = t - T[j-1]
        dist = distf + distb
        two = (distf/dist)*X[j] + (distb/dist)*X[j-1]
        v = (two - one)/(V*dT)
        pos = (two + one)/(2*V)
        poss[i] = pos
        A1[i], A2[i], A3[i] =  Act(F,dF,de,inve,D,pos,V,v)
        A1[i] *= dT
        A2[i] *= dT
        A3[i] *= dT
        one = two # move old posistion to new posistion
    end
    plot(poss)
    savefig("../Results/Graph2.png")
    println(sum(A1)/T[end])
    println(sum(A2)/T[end])
    println(sum(A3)/T[end])
    return(nothing)
end

@time main()
