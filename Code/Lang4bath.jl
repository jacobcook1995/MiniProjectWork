#!/usr/bin/env julia
# Lang4bath.jl
# A script to simulate euler-langevin trajectories about the steady state, and then
# calculate the entropy production
using Plots
using Roots
using DifferentialEquations
using SymEngine
import GR

# function to find the zeros of the function
function nullcline(F::Number,r::Number,f::Number,K::Number,Q::Number,k::Number,q::Number,kmin::Number,qmin::Number,high2low::Bool,Ne::Number)
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

# Vector of functions from MAP case
function fs(ps::AbstractVector)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin, A, B, S, W = symbols("K k Q q kmin qmin f r F Kmin Qmin A B S W")
    sym = Array{SymEngine.Basic}(4)
    sym[1] = k*S*r/(r+f*B^2) - K*A - kmin*A + Kmin*W
    sym[2] = q*S*r/(r+f*A^2) - Q*B - qmin*B + Qmin*W
    sym[3] = -k*S*r/(r + f*B^2) - q*S*r/(r + f*A^2) + kmin*A + qmin*B + F
    sym[4] = K*A + Q*B - Kmin*W - Qmin*W - F
    for i = 1:4
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        sym[i] = subs(sym[i], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
    end
    return sym
end

function f!(G::AbstractVector, x::AbstractVector, sym::Array{SymEngine.Basic,1})
    A, B, S, W = symbols("A B S W")
    for i = 1:4
        G[i] = subs(sym[i], A=>x[1], B=>x[2], S=>x[3], W=>x[4]) |> float
    end
    return G
end

# Diffusion matrix from MAP case
function Ds(ps::AbstractVector)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin, A, B, S, W = symbols("K k Q q kmin qmin f r F Kmin Qmin A B S W")
    e = Array{SymEngine.Basic,2}(4,4)
    e[1,2:4] = e[2,1] = e[2,3:4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    e[4,1] = -sqrt(K*A + Kmin*W)
    e[4,2] = -sqrt(Q*B + Qmin*W)
    e[4,4] = sqrt(F) # this is not as previously discussed, could be the error
    eT = transpose(e)
    D = e*eT
    D = inv(D)
    for j = 1:4
        for i = 1:4
           D[i,j] = subs(D[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
           D[i,j] = subs(D[i,j], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
       end
    end
    return D
end

function D!(D::AbstractArray, x::AbstractVector, par::Array{SymEngine.Basic,2})
    A, B, S, W = symbols("A B S W")
    # Make a symbolic version of the matrix, needs no input in this case
    for j = 1:4
        for i = 1:4
           D[i,j] = subs(par[i,j], A=>x[1], B=>x[2], S=>x[3], W=>x[4]) |> float
       end
    end
    return D
end

# Function to calculate the action of the discretised path in the MAP formulation
function EntProd(pathmin,tau,NM,ps)
    # probably easiest to calculate the entropy production at each point in the path
    F = ps[9]
    termsthi = zeros(4,4,NM)
    termsthif = zeros(4,4,NM)
    termsf = zeros(4,4,NM)
    Fqs = zeros(4,NM)
    Ffs = zeros(4,NM)
    Fs = zeros(NM)
    Acts = zeros(NM)
    Ents1 = zeros(NM)
    Ents2 = zeros(NM)
    h = [ 0.0; 0.0; 0.0; 0.0 ]
    thiv = [ 0.0; 0.0; 0.0; 0.0 ]
    d = [ 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 ]
    deltat = tau/NM
    # Remove fixed points from path
    path = pathmin[2:NM]
    fp = fs(ps)
    Dmin = Ds(ps)
    for i = 1:NM # This gives an extra contribution compared to the optimisation!
        if i == 1
            posA = (pathmin[1][1] + path[i][1])/2
            posB = (pathmin[1][2] + path[i][2])/2
            posS = (pathmin[1][3] + path[i][3])/2
            posW = (pathmin[1][4] + path[i][4])/2
        elseif i == NM
            posA = (path[i-1][1] + pathmin[NM+1][1])/2
            posB = (path[i-1][2] + pathmin[NM+1][2])/2
            posS = (path[i-1][3] + pathmin[NM+1][3])/2
            posW = (path[i-1][4] + pathmin[NM+1][4])/2
        else
            posA = (path[i-1][1] + path[i][1])/2
            posB = (path[i-1][2] + path[i][2])/2
            posS = (path[i-1][3] + path[i][3])/2
            posW = (path[i-1][4] + path[i][4])/2
        end
        h = f!(h, [posA; posB; posS; posW], fp)
        h[3] -= F
        h[4] += F
        d = D!(d, [posA; posB; posS; posW], Dmin)
        for j = 1:4
            if i == 1
                thiv[j] = (path[i][j] - pathmin[1][j])/deltat
            elseif i == NM
                thiv[j] = (pathmin[NM+1][j] - path[i-1][j])/deltat
            else
                thiv[j] = (path[i][j] - path[i-1][j])/(deltat)
            end
        end
        # Now calcualtion step
        for k = 1:4
            for j = 1:4
                termsthi[j,k,i] = thiv[j]*d[j,k]*thiv[k]
                termsthif[j,k,i] = h[j]*d[j,k]*thiv[k]
                termsf[j,k,i] = h[j]*d[j,k]*h[k]
            end
            Fqs[k,i] = thiv[k]*d[k,3]*F - thiv[k]*d[k,4]*F
            Ffs[k,i] = h[k]*d[k,3]*F - h[k]*d[k,4]*F
        end
        Fs[i] = F*d[3,3]*F + F*d[4,4]*F
    end
    # Now use these terms to calculate overall action along path and the entropy production
    for i = 1:NM
        Acts[i] = Fs[i]*(deltat/2)
        for k = 1:4
            for j = 1:4
                Acts[i] += termsthi[j,k,i]*(deltat/2)
                Acts[i] -= termsthif[j,k,i]*deltat
                Acts[i] += termsf[j,k,i]*(deltat/2)
                Ents1[i] += termsthif[j,k,i]*(2*deltat)
            end
            Acts[i] -= Fqs[k,i]*deltat
            Acts[i] += Ffs[k,i]*deltat
            Ents2[i] += Fqs[k,i]*(2*deltat) # comment this out to see the effect at somepoint
        end
    end

    return(Acts,Ents1,Ents2,termsthi,termsthif,termsf,Fqs,Ffs,Fs)
end

# f: function A = u[1], B = u[2], S = u[3]
function f1(du, u, p, t, K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin)
    du[1] = k*u[3]*r/(r + f*u[2]^2) - K*u[1] - kmin*u[1] + Kmin*u[4]
    du[2] = q*u[3]*r/(r + f*u[1]^2) - Q*u[2] - qmin*u[2] + Qmin*u[4]
    du[3] = F + kmin*u[1] + qmin*u[2] - k*u[3]*r/(r + f*u[2]^2) - q*u[3]*r/(r + f*u[1]^2)
    du[4] = K*u[1] + Q*u[2] - Kmin*u[4] - Qmin*u[4] - F
    return du
end

# g: noise function A = u[1], B = u[2], S = u[3]
function g1(du, u, p, t, K, k, Q, q, kmin, qmin, f, r, F, Ω, Kmin, Qmin)
    # Need to define a 4 by four matrix as I have 4 weiner processes and 4 dependant variables
    A = u[1]
    B = u[2]
    S = u[3]
    W = u[4]
    du[1,2:4] = du[2,1] = du[2,3:4] = du[3,4] = du[4,3] = 0
    du[1,1] = sqrt((k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)/Ω)
    du[2,2] = sqrt((q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)/Ω)
    du[3,1] = -sqrt((k*S*r/(r + f*B^2) + kmin*A)/Ω)
    du[3,2] = -sqrt((q*S*r/(r + f*A^2) + qmin*B)/Ω)
    du[3,3] = sqrt(F/Ω)
    du[4,1] = -sqrt((K*A + Kmin*W)/Ω)
    du[4,2] = -sqrt((Q*B + Qmin*W)/Ω)
    du[4,4] = sqrt(F/Ω)
    return du
end

function main()
    # General parameters
    Ω = 60 # system size, this is just a fudge to get my Euler-Maruyama algorithm (later) to work
    Ωr = Ω # fudge
    K = 1
    k = 1
    Q = 1
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1/((Ω/Ωr)^2) # Promoter switching
    r = 10
    F = 10*(Ω/Ωr)
    Kmin = 10.0^-20 # remains neligable though
    Qmin = 10.0^-20
    Ne = 150*(Ω/Ωr)
    T = 1
    high2low = true
    ps = [ K; k; Q; q; kmin; qmin; f; r; F; Kmin; Qmin ]
    n = 10000 # number of euler runs to do
    data = zeros(24,n) # preallocate memory
    # find start, mid and end points
    star, mid, fin = nullcline(F,r,f,K,Q,k,q,kmin,qmin,high2low,Ne)
    # make closed functions
    fc(du,u,p,t) = f1(du,u,p,t,K,k,Q,q,kmin,qmin,f,r,F,Kmin,Qmin)
    gc(du,u,p,t) = g1(du,u,p,t,K,k,Q,q,kmin,qmin,f,r,F,Ω,Kmin,Qmin)
    # First generate a run of the problem
    u₀ = star
    dt = (1/2)^(10)
    tspan = (0.0,T)
    prob = SDEProblem(fc, gc, u₀, tspan, noise_rate_prototype = zeros(4,4)) # SDEProblem
    Numb  = convert(Int64, floor(T/dt))
    for i = 1:n
        sol = DifferentialEquations.solve(prob, EM(), dt = dt) # To avoid namespace conflict with SymPy
        # Now would like to calculate both entropy production and the action along this path
        # now can use these paths to carry out a calculation of the action via MAP
        NM = size(sol[:],1)-1
        acts1, ents1, ents2, termsthi, termsthif, termsf, Fqs, Ffs, Fs = EntProd(sol,T,NM,ps)
        termsthi = squeeze(sum(termsthi,3),3)
        termsthi = squeeze(sum(termsthi,2),2)
        termsthif = squeeze(sum(termsthif,3),3)
        termsthif = squeeze(sum(termsthif,2),2)
        termsf = squeeze(sum(termsf,3),3)
        termsf = squeeze(sum(termsf,2),2)
        Fqs = squeeze(sum(Fqs,2),2)
        Ffs = squeeze(sum(Ffs,2),2)
        Fs = sum(Fs)
        data[1,i] = sum(acts1)/T
        data[2,i] = sum(ents1)/T
        data[3,i] = sum(ents2)/T
        data[4:7,i] = termsthi
        data[8:11,i] = termsthif
        data[12:15,i] = termsf
        data[16:19,i] = Fqs
        data[20:23,i] = Ffs
        data[24,i] = Fs
        if i % 100 == 0
            println(i)
        end
    end
    # Now need to write out this data to are file
    if length(ARGS) >= 1
        output_file = "../Results/$(ARGS[1]).csv"
        out_file = open(output_file, "w")
        # open file for writing
        for i = 1:size(data,2)
            line1 = "$(data[1,i]),$(data[2,i]),$(data[3,i]),$(data[4,i]),$(data[5,i]),$(data[6,i]),"
            line2 = "$(data[7,i]),$(data[8,i]),$(data[9,i]),$(data[10,i]),$(data[11,i]),$(data[12,i]),"
            line3 = "$(data[13,i]),$(data[14,i]),$(data[15,i]),$(data[16,i]),$(data[17,i]),$(data[18,i]),"
            line4 = "$(data[19,i]),$(data[20,i]),$(data[21,i]),$(data[22,i]),$(data[23,i]),$(data[24,i])\n"
            line = line1 * line2 * line3 * line4
            write(out_file, line)
        end
        # then close file
        close(out_file)
    end
end

@time main()
