#!/usr/bin/env julia
# Schlögltest.jl
using Roots

function main()
    k1 = 0.5 # k_{+1}[A]
    K1 = 3.0 # k_{-1}
    k2 = 1.0 # k_{+2}
    K2 = 1.0 # k_{-2}
    V = 10^8#0000#00
    B = 4.0
    f(x) = k1 - K1*x - k2*(x^3) + K2*B*(x^2)
    X = fzeros(f, 0, 10, no_pts = 1000)
    println(X)
    x = round(X[3]*V)
    println(x)
    R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
    println(R)
    S = 2*(((R[1]-R[2])^2) + ((R[3]-R[4])^2))/sum(R)
    println(S/V)
    println(27.758672288113193/(S/V))
    x = round(X[1]*V)
    println(x)
    R = [ k1*V, K1*x, k2*x*(x-1)*(x-2)/(V*V), K2*B*x*(x-1)/V]
    println(R)
    S = 2*(((R[1]-R[2])^2) + ((R[3]-R[4])^2))/sum(R)
    println(S/V)
    println(0.6732982218842689/(S/V))
    # Si = (2*(R[1]^2 + R[2]^2 + R[3]^2 + R[4]^2) + 4*(R[1]*R[4] + R[2]*R[3]))/sum(R)
    # println(Si/V)
    # Se = 4*(R[1]*R[2] + R[3]*R[4] + R[1]*R[3] + R[2]*R[4])/sum(R)
    # println(Se/V)
    return(nothing)
end

@time main()
