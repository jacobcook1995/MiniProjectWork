#!/usr/bin/env julia
# zeros_4bath.jl
# A script to find the zeros in terms of A so that the steady states and the
# inflex points can be found

# Putting the relevant imports in
using Optim
using Plots
using Roots
import GR # Need this to stop world age plotting error?

# Parameters
const Ω = 30 # system size, less important parameter now but should still retain
const ϕ = 0.1 # ratio ϕ = q/k
const k = 100 # steady state for A=k/K=1
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const K = k/Ω # K=k'
const Kmin = 10.0^-20
const q = k*ϕ # steady state for B=q/Q=1
const qmin = 10.0^-20
const Q = q/Ω # Q=q'
const Qmin = 10.0^-20
const f = 10 # Promoter switching
const r = 10
const F = 250 # removal rate
const N = 12000.0 # density within the system

# A = x[1], B = x[2], W = x[3], S = x[4]
function Wdot(x)
    wdot = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3] - F
    return(wdot)
end

function Zeros()
    g(x) = F./(1 .+ ((r .+ f*x.^2)./(ϕ*r .+ (f/ϕ)*((F/K .- x).^2)))) .+ K*x .- F
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
        Bs[i] = (1/ϕ)*((F/K) - As[i])
        Ss[i] = F/(k*r*(1/(r + f*Bs[i]^2) + ϕ/(r + f*As[i]^2)))
        Ws[i] = N - As[i] - Bs[i] - Ss[i]
    end
    print("$As,$Bs,$Ws,$Ss\n")
    vars = [ As Bs Ws Ss ]
    gr()
    A = 0:0.1:200
    gs = zeros(round(Int64,N+1))
    gs = g(A)
    plot(A,gs)
    savefig("../Results/GraphZeros.png")
    Wt = zeros(n)
    for i = 1:n
        Wt[i] = Wdot([ As[i] Bs[i] Ws[i] Ss[i] ])
    end
    print("$Wt\n")
    return(vars)
end

@time vars = Zeros()
