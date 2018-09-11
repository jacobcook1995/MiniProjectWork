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
function nullcline()
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
    # [1.12702, 8.87298, 15.0, 125.0]
    # [5.0, 5.0, 26.25, 113.75]
    # [8.87298, 1.12702, 15.0, 125.0]
    # now sub all parameters into dA
    Ω = 60 # system size, this is just a fudge to get my Euler-Maruyama algorithm (later) to work
    Ki = 1
    ki = 1
    Qi = 1
    qi = 11/15
    kmini = 0.5 # now reverse creation is an important process
    qmini = 0.1
    fi = 1000/((Ω/60)^2) # Promoter switching
    ri = 10000
    Fi = 10*(Ω/60)
    Kmini = 10.0^-3 # remains neligible though
    Qmini = 10.0^-3
    Ni = 150*(Ω/60)
    ps = [ Ki, ki, Qi, qi, kmini, Kmini, qmini, Qmini, fi, ri, Fi, Ni ]
    dA = subs(dA,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    dA = subs(dA,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # now find three stationary points in A
    guess = 0
    δ = 0.1
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
        if guess >= 10.0
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

@time nullcline()
