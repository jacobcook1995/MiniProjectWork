#!/usr/bin/env julia
# read4bath.jl
# A script to read in my file and find distinct points
using Plots
using Roots
using NLsolve
using SymEngine # Trying alternative substitution method to see if this is faster
ENV["GKSwstype"]="png" # ????????? this might stop the GKSterm error
import GR

# make a symbolic diffusion matrix
function Ds()
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,4,4)
    e[1,2:4] .= e[2,1] = e[2,3] = e[2,4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    e[4,1] = -sqrt(K*A + Kmin*W)
    e[4,2] = -sqrt(Q*B + Qmin*W)
    e[4,4] = sqrt(F) # this is not as previously discussed, could be the error
    # Now do the transformations required
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,4)
    b[1] = k*S*r/(r + f*B^2) - kmin*A - K*A + Kmin*W
    b[2] = q*S*r/(r + f*A^2) - qmin*B - Q*B + Qmin*W
    b[3] = -k*S*r/(r + f*B^2) - q*S*r/(r + f*A^2) + kmin*A + qmin*B + F
    b[4] = K*A + Q*B - (Kmin + Qmin)*W - F
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    the1, the2, the3, the4 = symbols("the1 the2 the3 the4")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the1*b[1] + the2*b[2] + the3*b[3] + the4*b[4]
    H += 0.5*the1*(D[1,1]*the1 + D[1,2]*the2 + D[1,3]*the3 + D[1,4]*the4)
    H += 0.5*the2*(D[2,1]*the1 + D[2,2]*the2 + D[2,3]*the3 + D[2,4]*the4)
    H += 0.5*the3*(D[3,1]*the1 + D[3,2]*the2 + D[3,3]*the3 + D[3,4]*the4)
    H += 0.5*the4*(D[4,1]*the1 + D[4,2]*the2 + D[4,3]*the3 + D[4,4]*the4)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    A, B, S, W = symbols("A B S W")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{SymEngine.Basic,1}(undef,4)
    Hx[1] = diff(H, A)
    Hx[2] = diff(H, B)
    Hx[3] = diff(H, S)
    Hx[4] = diff(H, W)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    the1, the2, the3, the4 = symbols("the1 the2 the3 the4")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{SymEngine.Basic,1}(undef,4)
    Hθ[1] = diff(H, the1)
    Hθ[2] = diff(H, the2)
    Hθ[3] = diff(H, the3)
    Hθ[4] = diff(H, the4)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    the1, the2, the3, the4 = symbols("the1 the2 the3 the4")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{SymEngine.Basic,2}(undef,4,4)
    Hθθ[1,1] = diff(Hθ[1],the1)
    Hθθ[1,2] = diff(Hθ[1],the2)
    Hθθ[1,3] = diff(Hθ[1],the3)
    Hθθ[1,4] = diff(Hθ[1],the4)
    Hθθ[2,1] = diff(Hθ[2],the1)
    Hθθ[2,2] = diff(Hθ[2],the2)
    Hθθ[2,3] = diff(Hθ[2],the3)
    Hθθ[2,4] = diff(Hθ[2],the4)
    Hθθ[3,1] = diff(Hθ[3],the1)
    Hθθ[3,2] = diff(Hθ[3],the2)
    Hθθ[3,3] = diff(Hθ[3],the3)
    Hθθ[3,4] = diff(Hθ[3],the4)
    Hθθ[4,1] = diff(Hθ[4],the1)
    Hθθ[4,2] = diff(Hθ[4],the2)
    Hθθ[4,3] = diff(Hθ[4],the3)
    Hθθ[4,4] = diff(Hθ[4],the4)
    return(Hθθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ follwed by X
function Hθxs()
    A, B, S, W = symbols("A B S W")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{SymEngine.Basic,2}(undef,4,4)
    Hθx[1,1] = diff(Hθ[1],A)
    Hθx[1,2] = diff(Hθ[1],B)
    Hθx[1,3] = diff(Hθ[1],S)
    Hθx[1,4] = diff(Hθ[1],W)
    Hθx[2,1] = diff(Hθ[2],A)
    Hθx[2,2] = diff(Hθ[2],B)
    Hθx[2,3] = diff(Hθ[2],S)
    Hθx[2,4] = diff(Hθ[2],W)
    Hθx[3,1] = diff(Hθ[3],A)
    Hθx[3,2] = diff(Hθ[3],B)
    Hθx[3,3] = diff(Hθ[3],S)
    Hθx[3,4] = diff(Hθ[3],W)
    Hθx[4,1] = diff(Hθ[4],A)
    Hθx[4,2] = diff(Hθ[4],B)
    Hθx[4,3] = diff(Hθ[4],S)
    Hθx[4,4] = diff(Hθ[4],W)
    return(Hθx)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y1, y2, y3, y4 = symbols("y1 y2 y3 y4")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1] + b[1]*Dmin[1,2]*b[2] + b[1]*Dmin[1,3]*b[3] + b[1]*Dmin[1,4]*b[4]
    num += b[2]*Dmin[2,1]*b[1] + b[2]*Dmin[2,2]*b[2] + b[2]*Dmin[2,3]*b[3] + b[2]*Dmin[2,4]*b[4]
    num += b[3]*Dmin[3,1]*b[1] + b[3]*Dmin[3,2]*b[2] + b[3]*Dmin[3,3]*b[3] + b[3]*Dmin[3,4]*b[4]
    num += b[4]*Dmin[4,1]*b[1] + b[4]*Dmin[4,2]*b[2] + b[4]*Dmin[4,3]*b[3] + b[4]*Dmin[4,4]*b[4]
    den = 0
    den += y1*Dmin[1,1]*y1 + y1*Dmin[1,2]*y2 + y1*Dmin[1,3]*y3 + y1*Dmin[1,4]*y4
    den += y2*Dmin[2,1]*y1 + y2*Dmin[2,2]*y2 + y2*Dmin[2,3]*y3 + y2*Dmin[2,4]*y4
    den += y3*Dmin[3,1]*y1 + y3*Dmin[3,2]*y2 + y3*Dmin[3,3]*y3 + y3*Dmin[3,4]*y4
    den += y4*Dmin[4,1]*y1 + y4*Dmin[4,2]*y2 + y4*Dmin[4,3]*y3 + y4*Dmin[4,4]*y4
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2, y3, y4 = symbols("y1 y2 y3 y4")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,4)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    c[3] = λ*y3 - b[3]
    c[4] = λ*y4 - b[4]
    ϑ = Array{SymEngine.Basic,1}(undef,4)
    ϑ[1] = Dmin[1,1]*c[1] + Dmin[1,2]*c[2] + Dmin[1,3]*c[3] + Dmin[1,4]*c[4]
    ϑ[2] = Dmin[2,1]*c[1] + Dmin[2,2]*c[2] + Dmin[2,3]*c[3] + Dmin[2,4]*c[4]
    ϑ[3] = Dmin[3,1]*c[1] + Dmin[3,2]*c[2] + Dmin[3,3]*c[3] + Dmin[3,4]*c[4]
    ϑ[4] = Dmin[4,1]*c[1] + Dmin[4,2]*c[2] + Dmin[4,3]*c[3] + Dmin[4,4]*c[4]
    return(ϑ)
end

# function to generate one of each symbolic object and to substitute parameter values into it
# this is a computationally intensive function, need to minimize the number of calls
function gensyms(ps::AbstractVector)
    # create symbolic objects
    ϑ = ϑs()
    λ = λs()
    Hθ = Hθs()
    Hθx = Hθxs()
    Hθθ = Hθθs()
    Hx = Hxs()
    H = Hs()
    # specify symbols that will be substituted for
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin = symbols("K k Q q kmin qmin f r F Kmin Qmin")
    # now perform substitutions
    λ = subs(λ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
    λ = subs(λ, qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
    H = subs(H, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
    H = subs(H, qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
    for i = 1:4
        ϑ[i] = subs(ϑ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        ϑ[i] = subs(ϑ[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
        Hθ[i] = subs(Hθ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        Hθ[i] = subs(Hθ[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
        Hx[i] = subs(Hx[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
        Hx[i] = subs(Hx[i], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
    end
    for j = 1:4
        for i = 1:4
            Hθx[i,j] = subs(Hθx[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
            Hθx[i,j] = subs(Hθx[i,j], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
            Hθθ[i,j] = subs(Hθθ[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5])
            Hθθ[i,j] = subs(Hθθ[i,j], qmin=>ps[6], f=>ps[7], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11])
        end
    end
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to generate the variables needed for a given algoritm iteration
function genvars(x::AbstractArray,λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int,Nmid::Int)
    # define neccesary symbols
    A, B, S, W, y1, y2, y3, y4 = symbols("A B S W y1 y2 y3 y4")
    # calculate velocities
    xprim = fill(NaN, NG+1, 4)
    for i = 2:NG
        for j = 1:4
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 4)
    ϑt = Array{SymEngine.Basic,1}(undef,4)
    for j = 1:4
        for i = 2:NG
            ϑt[j] = subs(ϑ[j], A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
            ϑs[i,j] = subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
        end
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# Function to calculate the action of a given path
function Ŝ(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int)
    S = zeros(NG-1)
    for j = 1:4
        S[1] += (3/(2*NG))*(xprim[2,j]*ϑs[2,j])
        # Not excluding the midpoint as the contribution is vanishing
        for i = 3:NG-1
            S[i-1] += (1/NG)*(xprim[i,j]*ϑs[i,j])
        end
        S[NG-1] += (3/(2*NG))*(xprim[NG,j]*ϑs[NG,j])
    end
    return(S)
end

# function to find the times of each point
function times(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::AbstractVector,x::AbstractArray,NM::Int)
    # Make discrete vector of time points
    t = zeros(NM+1)
    for i = 1:NM+1
        t[i] = (i-1)*ts[end]/NM
    end
    # Find index of first element greater than t[i] in ts
    inds = fill(0,NM+1)
    j = 1
    for i = 1:NM+1
        higher = false
        while higher == false
            if ts[j] >= t[i] || j == length(ts)
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they are fixed
    path = zeros(NM+1,4)
    path[1,:] = x[1,:]
    path[NM+1,:] = x[end,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        for j = 1:4
            x₀ = x[one,j]
            x₁ = x[two,j]
            path[i,j] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
        end
    end
    return(path)
end

# Vector of functions from MAP case
function f!(G, x, ps)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin = symbols("K k Q q kmin qmin f r F Kmin Qmin")
    sym = Array{SymEngine.Basic}(undef,4)
    sym[1] = k*x[3]*r/(r+f*x[2]^2) - K*x[1] - kmin*x[1] + Kmin*x[4]
    sym[2] = q*x[3]*r/(r+f*x[1]^2) - Q*x[2] - qmin*x[2] + Qmin*x[4]
    sym[3] = -k*x[3]*r/(r + f*x[2]^2) - q*x[3]*r/(r + f*x[1]^2) + kmin*x[1] + qmin*x[2] + F
    sym[4] = K*x[1] + Q*x[2] - Kmin*x[4] - Qmin*x[4] - F
    for i = 1:4
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        G[i] = subs(sym[i], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11]) |> float
    end
    return G
end
# Vector of functions from MAP case
function f1!(G, x, ps)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin = symbols("K k Q q kmin qmin f r F Kmin Qmin")
    sym = Array{SymEngine.Basic}(undef,4)
    sym[1] = k*x[3]*r/(r+f*x[2]^2) - kmin*x[1]
    sym[2] = q*x[3]*r/(r+f*x[1]^2) - qmin*x[2]
    sym[3] = F#-k*x[3]*r/(r + f*x[2]^2) - q*x[3]*r/(r + f*x[1]^2) + kmin*x[1] + qmin*x[2]
    sym[4] = K*x[1] + Q*x[2] - Kmin*x[4] - Qmin*x[4]
    for i = 1:4
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        G[i] = subs(sym[i], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11]) |> float
    end
    return G
end
# Vector of functions from MAP case
function f2!(G, x, ps)
    K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin = symbols("K k Q q kmin qmin f r F Kmin Qmin")
    sym = Array{SymEngine.Basic}(undef,4)
    sym[1] = K*x[1] - Kmin*x[4]
    sym[2] = Q*x[2] - Qmin*x[4]
    sym[3] = k*x[3]*r/(r + f*x[2]^2) + q*x[3]*r/(r + f*x[1]^2) - kmin*x[1] - qmin*x[2]#-F
    sym[4] = F
    for i = 1:4
        sym[i] = subs(sym[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
        G[i] = subs(sym[i], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11]) |> float
    end
    return G
end

# Diffusion matrix from MAP case
function D!(D, x, ps)
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,4,4)
    e[1,2:4] .= e[2,1] = e[2,3] = e[2,4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    e[4,1] = -sqrt(K*A + Kmin*W)
    e[4,2] = -sqrt(Q*B + Qmin*W)
    e[4,4] = sqrt(F) # this is not as previously discussed, could be the error
    # Now do the transformations required
    eT = transpose(e)
    D = e*eT
    D = inv(D)
    for j = 1:4
        for i = 1:4
           D[i,j] = subs(D[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], qmin=>ps[6], f=>ps[7])
           D[i,j] = subs(D[i,j], r=>ps[8], F=>ps[9], Kmin=>ps[10], Qmin=>ps[11], A=>x[1], B=>x[2], S=>x[3], W=>x[4]) |> float
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
    termsf11 = zeros(4,4,NM)
    termsf12 = zeros(4,4,NM)
    termsf22 = zeros(4,4,NM)
    Fqs = zeros(4,NM)
    Ffs = zeros(4,NM)
    Fs = zeros(NM)
    Acts = zeros(NM)
    Ents = zeros(NM)
    Ents2 = zeros(NM)
    Kins = zeros(NM)
    Pots = zeros(NM)
    Prods = zeros(NM)
    Flows = zeros(NM)
    h = [ 0.0; 0.0; 0.0; 0.0 ]
    h1 = [ 0.0; 0.0; 0.0; 0.0 ]
    h2 = [ 0.0; 0.0; 0.0; 0.0 ]
    thiv = [ 0.0; 0.0; 0.0; 0.0 ]
    d = [ 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 ]
    deltat = tau/NM
    # Remove fixed points from path
    path = pathmin[2:NM,:]
    for i = 1:NM # This gives an extra contribution compared to the optimisation!
        if i == 1
            posA = (pathmin[1,1] + path[i,1])/2
            posB = (pathmin[1,2] + path[i,2])/2
            posS = (pathmin[1,3] + path[i,3])/2
            posW = (pathmin[1,4] + path[i,4])/2
        elseif i == NM
            posA = (path[i-1,1] + pathmin[NM+1,1])/2
            posB = (path[i-1,2] + pathmin[NM+1,2])/2
            posS = (path[i-1,3] + pathmin[NM+1,3])/2
            posW = (path[i-1,4] + pathmin[NM+1,4])/2
        else
            posA = (path[i-1,1] + path[i,1])/2
            posB = (path[i-1,2] + path[i,2])/2
            posS = (path[i-1,3] + path[i,3])/2
            posW = (path[i-1,4] + path[i,4])/2
        end
        h = f!(h,[posA posB posS posW],ps)
        h[3] -= F
        h[4] += F
        d = D!(d,[posA posB posS posW],ps)
        h1 = f1!(h1,[posA posB posS posW],ps)
        h2 = f2!(h2,[posA posB posS posW],ps)
        for j = 1:4
            if i == 1
                thiv[j] = (path[i,j] - pathmin[1,j])/deltat
            elseif i == NM
                thiv[j] = (pathmin[NM+1,j] - path[i-1,j])/deltat
            else
                thiv[j] = (path[i,j] - path[i-1,j])/(deltat)
            end
        end
        # Now calculation step
        for k = 1:4
            for j = 1:4
                termsthi[j,k,i] = thiv[j]*d[j,k]*thiv[k]
                termsthif[j,k,i] = h[j]*d[j,k]*thiv[k]
                termsf[j,k,i] = h[j]*d[j,k]*h[k]
                termsf11[j,k,i] = h1[j]*d[j,k]*h1[k]
                termsf12[j,k,i] = h1[j]*d[j,k]*h2[k]
                termsf22[j,k,i] = h2[j]*d[j,k]*h2[k]
            end
            Fqs[k,i] = thiv[k]*d[k,3]*F - thiv[k]*d[k,4]*F
            Ffs[k,i] = h[k]*d[k,3]*F - h[k]*d[k,4]*F
        end
        Fs[i] = F*d[3,3]*F + F*d[4,4]*F
    end
    # Now use these terms to calculate overall action along path and the entropy production
    for i = 1:NM
        Acts[i] = Fs[i]*(deltat/2)
        Pots[i] = Fs[i]*(deltat/2)
        for k = 1:4
            for j = 1:4
                Acts[i] += termsthi[j,k,i]*(deltat/2)
                Kins[i] += termsthi[j,k,i]*(deltat/2)
                Acts[i] -= termsthif[j,k,i]*deltat
                Acts[i] += termsf[j,k,i]*(deltat/2)
                Pots[i] += termsf[j,k,i]*(deltat/2)
                Ents[i] += termsthif[j,k,i]*(2*deltat)
                Ents2[i] += termsthif[j,k,i]*(2*deltat)
                Flows[i] += termsf12[j,k,i]*(2*deltat)
                Prods[i] += termsf11[j,k,i]*deltat
                Prods[i] += termsf22[j,k,i]*deltat
            end
            Acts[i] -= Fqs[k,i]*deltat
            Acts[i] += Ffs[k,i]*deltat
            Ents[i] += Fqs[k,i]*(2*deltat)
            Ents2[i] -= Ffs[k,i]*(2*deltat)
            Pots[i] += Ffs[k,i]*deltat
        end
    end
    return(Acts,Ents,Kins,Pots,Prods,Flows,Ents2)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray,NG::Int,Nmid::Int)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s1 = zeros(Nmid)
    s2 = zeros(NG+2-Nmid)
    s1[1] = 0
    s2[1] = 0
    for i = 2:Nmid
        dA = x[i,1] - x[i-1,1]
        dB = x[i,2] - x[i-1,2]
        dS = x[i,3] - x[i-1,3]
        dW = x[i,4] - x[i-1,4]
        s1[i] = s1[i-1] + sqrt(dA^2 + dB^2 + dS^2 + dW^2) # Could probably drop the sqrts to speed up the code
    end
    for i = 2:(NG+2-Nmid)
        dA = x[i+Nmid-1,1] - x[i+Nmid-2,1]
        dB = x[i+Nmid-1,2] - x[i+Nmid-2,2]
        dS = x[i+Nmid-1,3] - x[i+Nmid-2,3]
        dW = x[i+Nmid-1,4] - x[i+Nmid-2,4]
        s2[i] = s2[i-1] + sqrt(dA^2 + dB^2 + dS^2 + dW^2) # Could probably drop the sqrts to speed up the code
    end
    # Divide total arc length into equal segments
    ls1 = zeros(Nmid)
    ls2 = zeros(NG+2-Nmid)
    for i = 1:Nmid
        ls1[i] = (i-1)*s1[end]/(Nmid-1)
    end
    for i = 1:(NG+2-Nmid)
        ls2[i] = (i-1)*s2[end]/(NG+1-Nmid)
    end
    # Find first index greater than a ls[i] for each i
    inds1 = fill(0,Nmid)
    j = 1
    for i = 1:Nmid
        higher = false
        while higher == false
            if s1[j] >= ls1[i] || j == Nmid
                inds1[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    inds2 = fill(0,NG+2-Nmid)
    j = 1
    for i = 1:(NG+2-Nmid)
        higher = false
        while higher == false
            if s2[j] >= ls2[i] || j == NG + 2 - Nmid
                inds2[i] = j + Nmid - 1
                higher = true
            else
                j += 1
            end
        end
    end
    # First do mid points and end points as they should be fixed
    disx = zeros(NG+1,4)
    disx[1,:] = x[1,:]
    disx[Nmid,:] = x[Nmid,:]
    disx[NG+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:Nmid-1
        one = inds1[i] - 1
        two = inds1[i]
        s₀ = s1[one]
        s₁ = s1[two]
        for j = 1:4
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls1[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end

    for i = Nmid+1:NG
        one = inds2[i+1-Nmid] - 1
        two = inds2[i+1-Nmid]
        s₀ = s2[one+1-Nmid]
        s₁ = s2[two+1-Nmid]
        for j = 1:4
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

function main()
    # create array to hold read in data
    points1 = Array{Float64}(undef,0,4)
    points2 = Array{Float64}(undef,0,4)
    # check if file is provided then read this data
    if length(ARGS) > 1
        println("Reading in $(ARGS[1])")
        open(ARGS[1], "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                comma = fill(0,3)
                j = 0
                L = length(line)
                for i = 1:L
                    if line[i] == ','
                        j += 1
                        comma[j] = i
                    end
                end
                A = parse(Float64, line[1:(comma[1] - 1)])
                B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
                S = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
                W = parse(Float64, line[(comma[3] + 1):L])
                points1 = vcat(points1, [ A B S W ])
            end
        end
        println("Reading in $(ARGS[2])")
        open(ARGS[2], "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                comma = fill(0,3)
                j = 0
                L = length(line)
                for i = 1:L
                    if line[i] == ','
                        j += 1
                        comma[j] = i
                    end
                end
                A = parse(Float64, line[1:(comma[1] - 1)])
                B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
                S = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
                W = parse(Float64, line[(comma[3] + 1):L])
                points2 = vcat(points2, [ A B S W ])
            end
        end
    else
        println("Need to provide 2 files to read")
        return(nothing)
    end
    # find NM for each data set
    NG1 = size(points1,1) - 1
    NG2 = size(points2,1) - 1
    if NG1 != NG2
        println("Files unequal length")
        return(nothing)
    end
    # think firstly I'd like to try and plot the path in 3D
    plot3d(points1[:,1], points1[:,2], points1[:,3], title = "Path 1", leg = false, camera = (10,70))
    savefig("../Results/Graph3D1.png")
    plot3d(points2[:,1], points2[:,2], points2[:,3], title = "Path 2", leg = false, camera = (10,70))
    savefig("../Results/Graph3D2.png")
    # make comnbined vector so both can be plotted
    points = zeros(NG1+1,4,2)
    for j = 1:4
        for i = 1:NG1+1
            points[i,j,1] = points1[i,j]
            points[i,j,2] = points2[i,j]
        end
    end
    plot3d(points[:,1,:], points[:,2,:], points[:,3,:], title = "Path 1&2", camera = (10,70))
    savefig("../Results/Graph3D.png")
    plot(points[:,4,:])
    savefig("../Results/Waste.png")

    # Now moving onto entropy production calculations
    Ω = 60 # system size
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/((Ω/60)^2) # Promoter switching
    r = 10.0
    F = 10.0*(Ω/60)
    Kmin = 10.0^-10
    Qmin = 10.0^-10

    Nmid = convert(Int64, ceil((NG1+1)/2))
    paras = [ K; k; Q; q; kmin; qmin; f; r; F; Kmin; Qmin ]
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(paras)
    x1, xprim1, λs1, ϑs1, λprim1 = genvars(points1,λ,ϑ,NG1,Nmid)
    x2, xprim2, λs2, ϑs2, λprim2 = genvars(points2,λ,ϑ,NG2,Nmid)
    plot(λs1)
    savefig("../Results/lambda1.png")
    plot(λs2)
    savefig("../Results/lambda2.png")
    λs1old = λs1[1:end]
    λs2old = λs2[1:end]
    S1 = Ŝ(x1,xprim1,λs1,ϑs1,λprim1,NG1)
    println("Action Path 1 = $(sum(S1))")
    S2 = Ŝ(x2,xprim2,λs2,ϑs2,λprim2,NG2)
    println("Action Path 2 = $(sum(S2))")
    # change lambdas to be appropriate for time discretisation
    star1 = star2 = false
    i = j = 0
    while star1 == false
        i += 1
        if λs1[i] > 10.0^-5
            star1 = true
        end
    end
    λs1[1:i-1] .= λs1[i]
    while star2 == false
        j += 1
        if λs2[j] > 10.0^-5
            star2 = true
        end
    end
    λs2[1:j-1] .= λs2[j]
    # second lot of loops to take care of the middle
    mid1 = mid2 = trig1 = trig2 = false
    ki = i
    m = j
    ks = ke = ms = me = 0
    while mid1 == false
        ki += 1
        if trig1 == false && λs1[ki] < 10.0^-5
            trig1 = true
            ks = ki
        elseif trig1 == true && λs1[ki] > 10.0^-5
            ke = ki
            mid1 = true
        end
    end
    λs1[ks:ke-1] .= (λs1[ks-1] + λs1[ke])/2
    while mid2 == false
        m += 1
        if trig2 == false && λs2[m] < 10.0^-5
            trig2 = true
            ms = m
        elseif trig2 == true && λs2[m] > 10.0^-5
            me = m
            mid2 = true
        end
    end
    λs2[ms:me-1] .= (λs2[ms-1] + λs2[me])/2
    λs1[end] = λs1[end-1]
    t1 = times(x1,xprim1,λs1,ϑs1,λprim1,NG1)
    λs2[end] = λs2[end-1]
    t2 = times(x2,xprim2,λs2,ϑs2,λprim2,NG2)
    plot(λs1-λs1old)
    savefig("../Results/lambdas1.png")
    plot(λs2-λs2old)
    savefig("../Results/lambdas2.png")
    plot(t1)
    savefig("../Results/times1.png")
    plot(t2)
    savefig("../Results/times2.png")
    # define midA, midB, midS, midW
    midA1 = x1[ks,1]
    midB1 = x1[ks,2]
    midS1 = x1[ks,3]
    midW1 = x1[ks,4]
    midA2 = x2[ms,1]
    midB2 = x2[ms,2]
    midS2 = x2[ms,3]
    midW2 = x2[ms,4]

    # Now rediscretise path into time discretisation
    NM1 = NM2 = 600
    path1 = timdis(t1,x1,NM1)
    # plot(path1[:,1],path1[:,2],path1[:,3])
    # savefig("../Results/NewPath1.png")
    path2 = timdis(t2,x2,NM2)
    # plot(path2[:,1],path2[:,2],path2[:,3])
    # savefig("../Results/NewPath2.png")
    # now can use these paths to carry out a calculation of the action via MAP
    path1r = path1[end:-1:1,:]
    path2r = path2[end:-1:1,:]
    ps = [ K; k; Q; q; kmin; qmin; f; r; F; Kmin; Qmin]
    acts1, ents1, kins1, pots1 = EntProd(path1,t1[end],NM1,ps)
    acts2, ents2, kins2, pots2 = EntProd(path2,t2[end],NM2,ps)
    acts1r, ents1r, kins1r, pots1r = EntProd(path1r,t1[end],NM1,ps)
    acts2r, ents2r, kins2r, pots2r = EntProd(path2r,t2[end],NM2,ps)
    points1 = [ 0.5*ents1, kins1, pots1, acts1]
    points2 = [ 0.5*ents2, kins2, pots2, acts2]
    println(sum(ents2)+sum(ents1r))
    println(sum(ents1)+sum(ents2r))
    println(sum(ents2)-sum(ents1))
    println(sum(ents1)-sum(ents2))
    println(sum(ents1))
    println(sum(ents2))
    println(sum(acts1))
    println(sum(acts2))
    # pone = plot(path1[1:300,1], points1, xaxis = "arc point", yaxis = "Action Contributions", marker = :auto, legend = false)
    # pone = scatter!(pone, [path1[1,1]], [0.0], seriescolor = :green)
    # pone = scatter!(pone, [midA1], [0.0], seriescolor = :orange)
    # pone = scatter!(pone, [path1[end,1]], [0.0], seriescolor = :red)
    # savefig("../Results/EntropyProduction1.png")
    # ptwo = plot(path2[1:300,1], points2, xaxis = "arc point", yaxis = "Action Contributions", marker = :auto, legend = false)
    # ptwo = scatter!(ptwo, [path2[1,1]], [0.0], seriescolor = :green)
    # ptwo = scatter!(ptwo, [midA2], [0.0], seriescolor = :orange)
    # ptwo = scatter!(ptwo, [path2[end,1]], [0.0], seriescolor = :red)
    # savefig("../Results/EntropyProduction2.png")
    # pone = plot(path1[:,1], path1[:,2], xaxis = "A", yaxis = "B", legend = false)
    # pone = scatter!(pone, [path1[1,1]], [path1[1,2]], seriescolor = :green)
    # pone = scatter!(pone, [midA1], [midB1], seriescolor = :orange)
    # pone = scatter!(pone, [path1[end,1]], [path1[end,2]], seriescolor = :red)
    # savefig("../Results/AvsB1.png")
    # pone = plot(path1[:,3], path1[:,4], xaxis = "S", yaxis = "W", legend = false)
    # pone = scatter!(pone, [path1[1,3]], [path1[1,4]], seriescolor = :green)
    # pone = scatter!(pone, [midS1], [midW1], seriescolor = :orange)
    # pone = scatter!(pone, [path1[end,3]], [path1[end,4]], seriescolor = :red)
    # savefig("../Results/SvsW1.png")
    # pone = plot(path2[:,1], path2[:,2], xaxis = "A", yaxis = "B", legend = false)
    # pone = scatter!(pone, [path2[1,1]], [path2[1,2]], seriescolor = :green)
    # pone = scatter!(pone, [midA2], [midB2], seriescolor = :orange)
    # pone = scatter!(pone, [path2[end,1]], [path2[end,2]], seriescolor = :red)
    # savefig("../Results/AvsB2.png")
    # pone = plot(path2[:,3], path2[:,4], xaxis = "S", yaxis = "W", legend = false)
    # pone = scatter!(pone, [path2[1,3]], [path2[1,4]], seriescolor = :green)
    # pone = scatter!(pone, [midS2], [midW2], seriescolor = :orange)
    # pone = scatter!(pone, [path2[end,3]], [path2[end,4]], seriescolor = :red)
    # savefig("../Results/SvsW2.png")
end

function shannon(point::Array{Float64,1},ps::Array{Float64,1})
    # should actually spell out what ps is here, so that
    # ps = [K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r]
    pa = ps[2]*ps[10]*point[3]/(ps[10]+ps[9]*(point[2])^2)
    pamin = ps[5]*point[1]
    pb = ps[4]*ps[10]*point[3]/(ps[10]+ps[9]*(point[1])^2)
    pbmin = ps[7]*point[2]
    da = ps[1]*point[1]
    damin = ps[6]*point[4]
    db = ps[3]*point[2]
    dbmin = ps[8]*point[4]
    F1 = (pa - pamin)
    F2 = (pb - pbmin)
    F3 = (da - damin)
    F4 = (db - dbmin)
    A1 = log(pa/pamin)
    A2 = log(pb/pbmin)
    A3 = log(da/damin)
    A4 = log(db/dbmin)
    S1 = F1*A1
    S2 = F2*A2
    S3 = F3*A3
    S4 = F4*A4
    S12 = S1 + S2
    S34 = S3 + S4
    S = S12 + S34
    return(S)
end

# function to do the same as graphs2 but for 4 species case
function graphs3()
    input_file1 = "../Results/1809/$(ARGS[1])1.csv"
    input_file2 = "../Results/1809/$(ARGS[1])2.csv"
    input_filep = "../Results/1809/$(ARGS[1])p.csv"
    points1 = Array{Float64,2}(undef,0,4)
    points2 = Array{Float64,2}(undef,0,4)
    ps = Array{Float64,1}(undef,0)
    # ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, Ne ]
    open(input_filep, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            p = parse(Float64, line)
            ps = vcat(ps, p)
        end
    end
    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = fill(0,3)
            L = length(line)
            i = 1
            for j = 1:L
                if line[j] == ','
                    comma[i] = j
                    i += 1
                end
            end
            A = parse(Float64, line[1:(comma[1] - 1)])
            B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
            S = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
            W = parse(Float64, line[(comma[3] + 1):L])
            points1 = vcat(points1, [ A B S W ])
        end
    end
    # now do second file
    open(input_file2, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = fill(0,3)
            L = length(line)
            i = 1
            for j = 1:L
                if line[j] == ','
                    comma[i] = j
                    i += 1
                end
            end
            A = parse(Float64, line[1:(comma[1] - 1)])
            B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
            S = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
            W = parse(Float64, line[(comma[3] + 1):L])
            points2 = vcat(points2, [ A B S W ])
        end
    end
    NG1 = size(points1,1) - 1
    NG2 = size(points2,1) - 1
    # right this is time discretised
    Nmid = convert(Int64, ceil((NG1+1)/2))
    # this is an inherently stupid step, to account for the different ps orderings
    # ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, Ne ]
    # ps2 = [K, k, Q, q, kmin, qmin, f, r, F, Kmin, Qmin]
    # this should be changed in future
    ps2 = [ ps[1], ps[2], ps[3], ps[4], ps[5], ps[7], ps[9], ps[10], ps[11], ps[6], ps[8] ]
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(ps2)
    x1, xprim1, λs1, ϑs1, λprim1 = genvars(points1,λ,ϑ,NG1,Nmid)
    x2, xprim2, λs2, ϑs2, λprim2 = genvars(points2,λ,ϑ,NG2,Nmid)
    λs1old = λs1[1:end]
    λs2old = λs2[1:end]
    S1 = Ŝ(x1,xprim1,λs1,ϑs1,λprim1,NG1)
    println("Action Path 1 = $(sum(S1))")
    S2 = Ŝ(x2,xprim2,λs2,ϑs2,λprim2,NG2)
    println("Action Path 2 = $(sum(S2))")
    # change lambdas to be appropriate for time discretisation
    star1 = star2 = false
    i = j = 0
    while star1 == false
        i += 1
        if λs1[i] > 10.0^-5
            star1 = true
        end
    end
    λs1[1:i-1] .= λs1[i]
    while star2 == false
        j += 1
        if λs2[j] > 10.0^-5
            star2 = true
        end
    end
    λs2[1:j-1] .= λs2[j]
    # second lot of loops to take care of the middle
    mid1 = mid2 = trig1 = trig2 = false
    ki = i
    m = j
    ks = ke = ms = me = 0
    while mid1 == false
        ki += 1
        if trig1 == false && λs1[ki] < 10.0^-5
            trig1 = true
            ks = ki
        elseif trig1 == true && λs1[ki] > 10.0^-5
            ke = ki
            mid1 = true
        end
    end
    λs1[ks:ke-1] .= (λs1[ks-1] + λs1[ke])/2
    while mid2 == false
        m += 1
        if trig2 == false && λs2[m] < 10.0^-5
            trig2 = true
            ms = m
        elseif trig2 == true && λs2[m] > 10.0^-5
            me = m
            mid2 = true
        end
    end
    λs2[ms:me-1] .= (λs2[ms-1] + λs2[me])/2
    λs1[end] = λs1[end-1]
    t1 = times(x1,xprim1,λs1,ϑs1,λprim1,NG1)
    λs2[end] = λs2[end-1]
    t2 = times(x2,xprim2,λs2,ϑs2,λprim2,NG2)
    println(t1[end])
    println(t2[end])
    # use this to time discretise the paths
    NM1 = NM2 = 600
    path1 = timdis(t1,points1,NM1) # high B to high A
    path2 = timdis(t2,points2,NM2) # high A to high B
    # now use these paths to find the relevant entropies
    acts1, ents1, kins1, pots1, prod1, flow1, entsF1 = EntProd(path1,t1[end],NM1,ps2)
    acts2, ents2, kins2, pots2, prod2, flow2, entsF2 = EntProd(path2,t2[end],NM2,ps2)
    println(sum(acts1))
    println(sum(acts2))
    # find segment centers to plot against
    segcent1 = zeros(NM1)
    segcent2 = zeros(NM2)
    for i = 1:length(segcent1)
        segcent1[i] = (path1[i,1] + path1[i+1,1])/2
    end
    for i = 1:length(segcent2)
        segcent2[i] = (path2[i,1] + path2[i+1,1])/2
    end
    S1 = shannon(points1[1,:],ps) # high B
    S2 = shannon(points2[1,:],ps) # high A
    # factor 2 exists to account for difference between entropy and entropic contribution to action
    plot(segcent1,2*NM1*ents1/t1[end],title="B => A",label="Entropy Production",xlabel="A",ylabel="\\Delta S")
    # plot!(segcent1,2*NM1*flow1/t1[end],label="''Entropy Flow''")
    # plot!(segcent1,2*NM1*prod1/t1[end],label="''Entropy Production''")
    # plot!(segcent1,2*NM1*(prod1-flow1)/t1[end],label="Production-Flow")
    plot!(segcent1,2*NM1*entsF1/t1[end],label="Ent Prod Reverse F")
    scatter!([segcent1[1]], [0.0], seriescolor = :green, label="")
    scatter!([segcent1[end]], [0.0], seriescolor = :red, label="")
    savefig("../Results/ents1.png")
    hline!([S1],label="high B Schnakenberg")
    hline!([S2],label="high A Schnakenberg")
    savefig("../Results/ents1sch.png")
    plot(segcent2,2*NM2*ents2/t2[end],title="A => B",label="Entropy Production",xlabel="A",ylabel="\\Delta S")
    # plot!(segcent2,2*NM2*flow2/t2[end],label="''Entropy Flow''")
    # plot!(segcent2,2*NM2*prod2/t2[end],label="''Entropy Production''")
    # plot!(segcent2,2*NM2*(prod2-flow2)/t2[end],label="Production-Flow")
    plot!(segcent2,2*NM2*entsF2/t2[end],label="Ent Prod Reverse F")
    scatter!([segcent2[1]], [0.0], seriescolor = :green, label="")
    scatter!([segcent2[end]], [0.0], seriescolor = :red, label="")
    savefig("../Results/ents2.png")
    hline!([S1],label="high B Schnakenberg")
    hline!([S2],label="high A Schnakenberg",legend=false)
    savefig("../Results/ents2sch.png")
    return(nothing)
end

# @time main()
@time graphs3()
