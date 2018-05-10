#!/usr/bin/env julia
# Lang3bath.jl
# A script to simulate euler-langevin trajectories about the steady state, and then
# calculate the entropy production
using Plots
using Roots

# function to find the zeros of the function
function nullcline(F::Number,r::Number,f::Number,K::Number,Q::Number,k::Number,q::Number,kmin::Number,qmin::Number,high2low::Bool)
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
    for i = 1:n
        Bs[i] = (F - K*As[i])/Q
        Ss[i] = (1/(k*r))*(r + f*((F - K*As[i])/Q)^2)*(K + kmin)*As[i]
    end
    sad = [ As[2]; Bs[2]; Ss[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ss[1] ]
        ss2 = [ As[3]; Bs[3]; Ss[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ss[3] ]
        ss2 = [ As[1]; Bs[1]; Ss[1] ]
    end
    print("$(ss1)\n")
    print("$(sad)\n")
    print("$(ss2)\n")
    return (ss1,sad,ss2)
end

function main()
    # General parameters
    K = 1
    k = 1
    Q = 1
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1 # Promoter switching
    r = 10
    F = 10
    high2low = true # somewhat meaningless in this script

    # find start, mid and end points
    star, mid, fin = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low)
end

@time main()
