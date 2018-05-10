#!/usr/bin/env julia
# Lang3bath.jl
# A script to simulate euler-langevin trajectories about the steady state, and then
# calculate the entropy production
using Plots
using Roots
using DifferentialEquations
using SymEngine
import GR

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

# Vector of functions from MAP case
function f!(G, x, ps)
    K, k, Q, q, kmin, qmin, f, r, F = symbols("K k Q q kmin qmin f r F")
    sym = Array{SymEngine.Basic}(3)
    sym[1] = k*x[3]*r/(r+f*x[2]^2) - K*x[1] - kmin*x[1]
    sym[2] = q*x[3]*r/(r+f*x[1]^2) - Q*x[2] - qmin*x[2]
    sym[3] = -k*x[3]*r/(r + f*x[2]^2) - q*x[3]*r/(r + f*x[1]^2) + kmin*x[1] + qmin*x[2] + F
    for i = 1:3
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        G[i] = subs(sym[i], r=>ps[8], F=>ps[9]) |> float
    end
    return G
end

# Diffusion matrix from MAP case
function D!(D, x, ps)
    A, B, S, K, k, Q, q, kmin, qmin, f, r, F = symbols("A B S K k Q q kmin qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(3,3)
    e[1,2:3] = e[2,1] = e[2,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    eT = transpose(e)
    D = e*eT
    D = inv(D)
    for j = 1:3
        for i = 1:3
           D[i,j] = subs(D[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
           D[i,j] = subs(D[i,j], r=>ps[8], F=>ps[9], A=>x[1], B=>x[2], S=>x[3]) |> float
       end
    end
    return D
end

# Function to calculate the action of the discretised path in the MAP formulation
function EntProd(pathmin,tau,NM,ps)
    # probably easiest to calculate the entropy production at each point in the path
    F = ps[9]
    termsthi = zeros(3,3,NM)
    termsthif = zeros(3,3,NM)
    termsf = zeros(3,3,NM)
    Fqs = zeros(3,NM)
    Ffs = zeros(3,NM)
    Fs = zeros(NM)
    Acts = zeros(NM)
    Ents = zeros(NM)
    h = [ 0.0; 0.0; 0.0 ]
    thiv = [ 0.0; 0.0; 0.0 ]
    d = [ 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0 ]
    deltat = tau/NM
    # Remove fixed points from path
    path = pathmin[2:NM]
    for i = 1:NM # This gives an extra contribution compared to the optimisation!
        if i == 1
            posA = (pathmin[1][1] + path[i][1])/2
            posB = (pathmin[1][2] + path[i][2])/2
            posS = (pathmin[1][3] + path[i][3])/2
        elseif i == NM
            posA = (path[i-1][1] + pathmin[NM+1][1])/2
            posB = (path[i-1][2] + pathmin[NM+1][2])/2
            posS = (path[i-1][3] + pathmin[NM+1][3])/2
        else
            posA = (path[i-1][1] + path[i][1])/2
            posB = (path[i-1][2] + path[i][2])/2
            posS = (path[i-1][3] + path[i][3])/2
        end
        h = f!(h, [posA posB posS], ps)
        h[3] -= F
        d = D!(d, [posA posB posS], ps)
        for j = 1:3
            if i == 1
                thiv[j] = (path[i][j] - pathmin[1][j])/deltat
            elseif i == NM
                thiv[j] = (pathmin[NM+1][j] - path[i-1][j])/deltat
            else
                thiv[j] = (path[i][j] - path[i-1][j])/(deltat)
            end
        end
        # Now calcualtion step
        for k = 1:3
            for j = 1:3
                termsthi[j,k,i] = thiv[j]*d[j,k]*thiv[k]
                termsthif[j,k,i] = h[j]*d[j,k]*thiv[k]
                termsf[j,k,i] = h[j]*d[j,k]*h[k]
            end
            Fqs[k,i] = thiv[k]*d[k,3]*F
            Ffs[k,i] = h[k]*d[k,3]*F
        end
        Fs[i] = F*d[3,3]*F
    end
    # Now use these terms to calculate overall action along path and the entropy production
    for i = 1:NM
        Acts[i] = Fs[i]*(deltat/2)
        for k = 1:3
            for j = 1:3
                Acts[i] += termsthi[j,k,i]*(deltat/2)
                Acts[i] -= termsthif[j,k,i]*deltat
                Acts[i] += termsf[j,k,i]*(deltat/2)
                Ents[i] += termsthif[j,k,i]*(2*deltat)
            end
            Acts[i] -= Fqs[k,i]*deltat
            Acts[i] += Ffs[k,i]*deltat
            Ents[i] += Fqs[k,i]*(2*deltat) # comment this out to see the effect at somepoint
        end
    end
    return(Acts,Ents)
end

# f: function A = u[1], B = u[2], S = u[3]
function f1(du, u, p, t, K, k, Q, q, kmin, qmin, f, r, F)
    du[1] = k*u[3]*r/(r + f*u[2]^2) - K*u[1] - kmin*u[1]
    du[2] = q*u[3]*r/(r + f*u[1]^2) - Q*u[2] - qmin*u[2]
    du[3] = F + kmin*u[1] + qmin*u[2] - k*u[3]*r/(r + f*u[2]^2) - q*u[3]*r/(r + f*u[1]^2)
    return du
end

# g: noise function A = u[1], B = u[2], S = u[3]
function g1(du, u, p, t, K, k, Q, q, kmin, qmin, f, r, F, Ω)
    # Need to define a 4 by four matrix as I have 4 weiner processes and 4 dependant variables
    A = u[1]
    B = u[2]
    S = u[3]
    du[1,2:3] = du[2,1] = du[2,3] = 0
    du[1,1] = sqrt((k*S*r/(r + f*B^2) + kmin*A + K*A)/Ω)
    du[2,2] = sqrt((q*S*r/(r + f*A^2) + qmin*B + Q*B)/Ω)
    du[3,1] = -sqrt((k*S*r/(r + f*B^2) + kmin*A)/Ω)
    du[3,2] = -sqrt((q*S*r/(r + f*A^2) + qmin*B)/Ω)
    du[3,3] = sqrt(F/Ω)
    return du
end

function main()
    # General parameters
    Ω = 60 # system size
    K = 1
    k = 1
    Q = 1
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1/((Ω/60)^2) # Promoter switching
    r = 10
    F = 10*(Ω/60)
    high2low = true # somewhat meaningless in this script
    T = 15#150 # time period studied

    # find start, mid and end points
    star, mid, fin = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low)
    # make closed functions
    fc(du,u,p,t) = f1(du,u,p,t,K,k,Q,q,kmin,qmin,f,r,F)
    gc(du,u,p,t) = g1(du,u,p,t,K,k,Q,q,kmin,qmin,f,r,F,Ω)
    # First generate a run of the problem
    u₀ = star
    dt = (1/2)^(10)
    tspan = (0.0,T)
    prob = SDEProblem(fc, gc, u₀, tspan, noise_rate_prototype = zeros(3,3)) # SDEProblem
    sol = DifferentialEquations.solve(prob, EM(), dt = dt) # To avoid namespace conflict with SymPy
    plot(sol)
    savefig("../Results/ProbablyNotAGoodGraph.png")
    # Now would like to calculate both entropy production and the action along this path
    # now can use these paths to carry out a calculation of the action via MAP
    ps = [ K; k; Q; q; kmin; qmin; f; r; F ]
    NM = size(sol[:],1)-1
    acts1, ents1 = EntProd(sol,T,NM,ps)
    plot(acts1)
    savefig("../Results/PossiblyMoreUseful.png")
    plot(ents1)
    savefig("../Results/PossiblyMoreUseful2.png")
end

@time main()
