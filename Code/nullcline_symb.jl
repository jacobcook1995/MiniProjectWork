#!/usr/bin/env julia
# nullcline_symb.jl
# A script to find the symbolic equations for the nullclines and the steady states
using SymPy
using Plots
import GR # Need this to stop world age plotting error?

# my function
function null()
    # Define my neccesary symbols
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N = symbols("A,B,S,W,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,N")
    # nullcline for W in terms of A and B
    Ss = N - A - B - W
    # write out expressions
    e1 = k*S*r/(r + f*B^2) - K*A - kmin*A + Kmin*W
    e2 = q*S*r/(r + f*A^2) - Q*B - qmin*B + Qmin*W
    e3 = -k*S*r/(r + f*B^2) - q*S*r/(r + f*A^2) + kmin*A + qmin*B + F
    e4 = K*A + Q*B - Kmin*W - Qmin*W - F
    # sub in new value for S
    e1 = subs(e1, S=>Ss)
    e2 = subs(e2, S=>Ss)
    e3 = subs(e3, S=>Ss)
    e4 = subs(e4, S=>Ss)
    # solve e4 for W
    Ws = solve(e4, W)[1] #[1] converts from symbolic array to symbolic expression
    # sub in new value for W
    e1 = subs(e1, W=>Ws)
    e2 = subs(e2, W=>Ws)
    e3 = subs(e3, W=>Ws)
    # solve e2 for B
    Bs = solve(e2, B)[1]
    # sub in new value for B
    e1 = subs(e1, B=>Bs, k => 15, kmin =>5, q=>11, qmin=>1, K=>1.5, Kmin=>0.5, Q=>1.1, Qmin=>0.1, F=>250, r=>1000, f=>20, N=>1000)
    e3 = subs(e3, B=>Bs, k => 15, kmin =>5, q=>11, qmin=>1, K=>1.5, Kmin=>0.5, Q=>1.1, Qmin=>0.1, F=>250, r=>1000, f=>20, N=>1000)
    points = zeros(10000,2)
    for i = 1:10000
        points[i,1] = subs(e1, A=>(i/10))
        points[i,2] = subs(e1, A=>(i/10))
    end
    plot(points)
    savefig("../Results/GrApH.png")
end

@time null()
