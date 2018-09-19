#!/usr/bin/env julia
# stead4bath.jl
# A script that uses machine algebra to find the steady state of the 4 species
# bistable toggle switch
#
# Author: Jacob Cook
# Date: September 2018

using SymPy
using Plots
using Roots
import GR

# function to find nullclines
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N ]
# maxr = range to look over, δ step size in this range
function nullcline(ps::Array{Float64,1},maxr::Float64,δ::Float64)
    # define symbols
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N = symbols("A,B,S,W,K k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,N")
    # equation for S
    Se = N - W - A - B
    # equation for W
    We = (K*A + Q*B - F)/(Qmin + Kmin)
    # eliminate W from equation for S
    Se = subs(Se,W=>We)
    # equation for dB/dt
    dB = (q*S*r)/(r + f*A^2) - qmin*B - Q*B + Qmin*W
    # sub in S and W expressions
    dB = subs(dB,S=>Se,W=>We)
    # solve this expression for B
    Bear = solve(dB,B)
    # remove from being in array form
    Be = Bear[1]
    # now sub this into Se and We
    Se = subs(Se,B=>Be)
    We = subs(We,B=>Be)
    # equation for dA/dt
    dA = (k*S*r)/(r + f*B^2) - kmin*A - K*A + Kmin*W
    # sub in B, S and W expressions
    dA = subs(dA,B=>Be,S=>Se,W=>We)
    # then sub in parameters
    dA = subs(dA,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    dA = subs(dA,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # now find three stationary points in A
    guess = 0
    three = false
    As = Array{Float64,1}(undef,0)
    n = 0
    while three == false
        a = convert(Float64,nsolve(dA,guess))
        # add to vector if not already found
        if (a in As) == false
            As = vcat(As,a)
            n += 1
        end
        if n == 3
            three = true
        end
        guess += δ
        if guess >= maxr
            println("Could not find three stationary points in range (0,$(guess)) with step $(δ)")
            println(As)
            error()
        end
    end
    # sub parameters into Be, Se, We
    Be = subs(Be,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Be = subs(Be,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    Se = subs(Se,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Se = subs(Se,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    We = subs(We,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    We = subs(We,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # then sub the various As values to get final values
    Bs = zeros(3)
    Ss = zeros(3)
    Ws = zeros(3)
    for i = 1:3
        Bs[i] = subs(Be,A=>As[i]) |> float
        Ss[i] = subs(Se,A=>As[i]) |> float
        Ws[i] = subs(We,A=>As[i]) |> float
    end
    # finally put into states
    ss1 = [As[1],Bs[1],Ss[1],Ws[1]]
    sad = [As[2],Bs[2],Ss[2],Ws[2]]
    ss2 = [As[3],Bs[3],Ss[3],Ws[3]]
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

function main()
    # General parameters
    Ω = 60 # system size, this is just a fudge to get my Euler-Maruyama algorithm (later) to work
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11.0/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/((Ω/60)^2) # Promoter switching
    r = 10.0
    F = 10.0*(Ω/60)
    Kmin = 10.0^-10 # remains neligable though
    Qmin = 10.0^-10
    N = 150.0*(Ω/60)
    ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N ]
    nullcline(ps,10.0,0.1)
    return(nothing)
end

@time main()
