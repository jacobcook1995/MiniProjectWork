#!/usr/bin/env julia
# ODEent.jl
# A script to calculate the entropy production from the ODE model at steady state

using Roots

# function to find the zeros of the function
function nullcline(F::Float64,r::Float64,f::Float64,K::Float64,Q::Float64,k::Float64,
                    q::Float64,kmin::Float64,qmin::Float64,high2low::Bool,Ne::Float64)
    g(x) = (K + kmin)*(q/k)*((r + f*((F - K*x)/Q)^2)/(r + f*x^2))*x - (qmin + Q)*(F - K*x)/Q
    three = false
    n = 0
    As = []
    while three == false
        As = fzeros(g, 0, 2*F/K, order = 1)
        n = length(As)
        if n == 3
            three = true
        end
    end
    Bs = zeros(n)
    Ss = zeros(n)
    Ws = zeros(n)
    for i = 1:n
        Bs[i] = (F - K*As[i])/Q
        Ss[i] = (1/(k*r))*(r + f*((F - K*As[i])/Q)^2)*(K + kmin)*As[i]
        Ws[i] = Ne - As[i] - Bs[i] - Ss[i]
    end
    sad = [ As[2]; Bs[2]; Ss[2]; Ws[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ss[1]; Ws[1] ]
        ss2 = [ As[3]; Bs[3]; Ss[3]; Ws[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ss[3]; Ws[3] ]
        ss2 = [ As[1]; Bs[1]; Ss[1]; Ws[1] ]
    end
    print("$(ss1)\n")
    print("$(sad)\n")
    print("$(ss2)\n")
    return (ss1,sad,ss2)
end

# function to calculate fluxes along both paths
function fluxes(point::Array{Float64,1},F::Float64,r::Float64,f::Float64,K::Float64,
                Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    F11 = K*point[1] - Kmin*point[4]
    F12 = k*point[3]*r/(r + f*point[2]^2) - kmin*point[1]
    F21 = Q*point[2] - Qmin*point[4]
    F22 = q*point[3]*r/(r + f*point[1]^2) - qmin*point[2]
    F1 = F11
    F2 = F21
    return(F1,F2)
end

# function to take in fluxes and calculate the steady state entropy production
function ent(point::Array{Float64,1},F1::Float64,F2::Float64,r::Float64,f::Float64,K::Float64,
                Q::Float64,k::Float64,q::Float64,kmin::Float64,qmin::Float64,Kmin::Float64,Qmin::Float64)
    S = 0
    S += F1*log(K*point[1]/(Kmin*point[4]))
    S += F1*log(k*point[3]*r/((r + f*point[2]^2)*kmin*point[1]))
    S += F2*log(Q*point[2]/(Qmin*point[4]))
    S += F2*log(q*point[3]*r/((r + f*point[1]^2)*qmin*point[2]))
    return(S)
end

function main()
    # General parameters
    Ω = 1 # Ω = 1, gMAP parameterisation
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/(Ω^2) # Promoter switching
    r = 10.0
    F = 10.0*Ω
    Kmin = 10.0^-20 # remains neligable though
    Qmin = 10.0^-20
    Ne = 150.0*Ω # number of elements in the system
    high2low = true
    star, mid, fin = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low,Ne)
    F1s, F2s = fluxes(star,F,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    Ss = ent(star,F1s,F2s,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println(Ss)
    F1m, F2m = fluxes(mid,F,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    Sm = ent(mid,F1m,F2m,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println(Sm)
    F1f, F2f = fluxes(fin,F,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    Sf = ent(fin,F1f,F2f,r,f,K,Q,k,q,kmin,qmin,Kmin,Qmin)
    println(Sf)

end

@time main()
