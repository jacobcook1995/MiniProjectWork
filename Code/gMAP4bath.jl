#!/usr/bin/env julia
# gMAP4bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to this
# time done for the 4 species case

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
#using SymPy
using SymEngine # Trying alternative substitution method to see if this is faster
import GR # Need this to stop world age plotting error?

# first should make generic functions that take in a full set of parameters, that
# can perform the machine alegbra and then sub in the parameter set to generate
# numeric values for each of the instances considered

# make a symbolic diffusion matrix
function Ds()
    #A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A,B,S,W,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F")
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{Sym}(4,4)
    e[1,2:4] = e[2,1] = e[2,3:4] = e[3,4] = e[4,3] = 0
    e[1,1] = sqrt(k*S*r/(r + f*B^2) + kmin*A + K*A + Kmin*W)
    e[2,2] = sqrt(q*S*r/(r + f*A^2) + qmin*B + Q*B + Qmin*W)
    e[3,1] = -sqrt(k*S*r/(r + f*B^2) + kmin*A)
    e[3,2] = -sqrt(q*S*r/(r + f*A^2) + qmin*B)
    e[3,3] = sqrt(F)
    e[4,1] = -sqrt(K*A + Kmin*W)
    e[4,2] = -sqrt(Q*B + Qmin*W)
    e[4,4] = sqrt(F)
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
    #A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A,B,S,W,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F")
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{Sym}(4)
    b[1] = k*S*r/(r + f*B^2) - kmin*A - K*A + Kmin*W
    b[2] = q*S*r/(r + f*A^2) - qmin*B - Q*B + Qmin*W
    b[3] = -k*S*r/(r + f*B^2) - q*S*r/(r + f*A^2) + kmin*A + qmin*B + F
    b[4] = K*A + Q*B - (Kmin + Qmin)*W - F
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    #the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
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
    #A, B, S, W = symbols("A,B,S,W")
    A, B, S, W = symbols("A B S W")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{Sym}(4)
    Hx[1] = diff(H, A)
    Hx[2] = diff(H, B)
    Hx[3] = diff(H, S)
    Hx[4] = diff(H, W)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    #the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
    the1, the2, the3, the4 = symbols("the1 the2 the3 the4")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{Sym}(4)
    Hθ[1] = diff(H, the1)
    Hθ[2] = diff(H, the2)
    Hθ[3] = diff(H, the3)
    Hθ[4] = diff(H, the4)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    #the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
    the1, the2, the3, the4 = symbols("the1 the2 the3 the4")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{Sym}(4,4)
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
    # A, B, S, W = symbols("A,B,S,W")
    A, B, S, W = symbols("A B S W")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{Sym}(4,4)
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
    #y1, y2, y3, y4 = symbols("y1,y2,y3,y4")
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
    #y1, y2, y3, y4 = symbols("y1,y2,y3,y4")
    y1, y2, y3, y4 = symbols("y1 y2 y3 y4")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{Sym}(4)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    c[3] = λ*y3 - b[3]
    c[4] = λ*y4 - b[4]
    ϑ = Array{Sym}(4)
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
    # specify symbols that will be substituted for
    #K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F")
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = symbols("K k Q q kmin Kmin qmin Qmin f r F")
    # now perform substitutions
    ϑ = subs(ϑ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    ϑ = subs(ϑ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    λ = subs(λ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    λ = subs(λ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    Hθ = subs(Hθ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    Hθ = subs(Hθ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    Hθx = subs(Hθx, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    Hθx = subs(Hθx, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    Hθθ = subs(Hθθ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    Hθθ = subs(Hθθ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    Hx = subs(Hx, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6]) |> Sym
    Hx = subs(Hx, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11]) |> Sym
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx)
end

# function to find the zeros of the function
function nullcline(F::Number,r::Number,f::Number,ϕ::Number,K::Number,k::Number,Ne::Int,high2low::Bool)
    # Define read in parameters
    g(x) = F./(1 + ((r + f*x.^2)./(ϕ*r + (f/ϕ)*((F/K - x).^2)))) + K.*x - F
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
    print(ss1)
    print("\n")
    print(sad)
    print("\n")
    print(ss2)
    print("\n")
    return (ss1,sad,ss2)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray,NG::Int)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s = zeros(NG+1)
    s[1] = 0
    for i = 2:NG+1
        dA = x[i,1] - x[i-1,1]
        dB = x[i,2] - x[i-1,2]
        dS = x[i,3] - x[i-1,3]
        dW = x[i,4] - x[i-1,4]
        s[i] = s[i-1] + sqrt(dA^2 + dB^2 + dS^2 + dW^2)
    end
    # Divide total arc length into equal segments
    ls = zeros(NG+1)
    for i = 1:NG+1
        ls[i] = (i-1)*s[end]/(NG)
    end
    # Find first index greater than a ls[i] for each i
    inds = fill(0,NG+1)
    j = 1
    for i = 1:NG+1
        higher = false
        while higher == false
            if s[j] >= ls[i] || j == NG+1
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they should be fixed
    disx = zeros(NG+1,4)
    disx[1,:] = x[1,:]
    disx[NG+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NG
        one = inds[i] - 1
        two = inds[i]
        s₀ = s[one]
        s₁ = s[two]
        for j = 1:4
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

# function to generate the variables needed for a given algoritm iteration
function genvars(x::AbstractArray, λ::SymPy.Sym, ϑ::SymPy.Sym, NG::Int)
    # define neccesary symbols
    #A, B, S, W, y1, y2, y3, y4 = symbols("A,B,S,W,y1,y2,y3,y4")
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
    for i = 2:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4]) |> Sym
        λt = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> Sym
        λs[i] = N(λt)
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 4)
    for i = 2:NG
        ϑt = ϑ # temporary ϑ to avoid altering master one
        ϑt = subs(ϑt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4]) |> Sym
        ϑs[i,:] = subs(ϑt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# function to be solved by NLsolve
function g!(F::AbstractArray, x::AbstractArray, C::AbstractVector, K::AbstractArray, xi::AbstractArray, NG::Int)
    # Start point
    F[1,:] = x[1,:] - xi[1,:]
    # first path
    for i = 2:NG
        for j = 1:4
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # end point
    F[end,:] = x[end,:] - xi[2,:]
    return(F)
end

# function to solve the system of linear equations
function linsys(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,
                Hx::SymPy.Sym,Hθ::SymPy.Sym,Hθθ::SymPy.Sym,Hθx::SymPy.Sym,Δτ::Number,NG::Int)
    # define relevant symbols
    #A, B, S, W, the1, the2, the3, the4 = symbols("A,B,S,W,the1,the2,the3,the4")
    A, B, S, W, the1, the2, the3, the4 = symbols("A B S W the1 the2 the3 the4")
    # Make array to store fixed points
    xi = fill(NaN, 2, 4)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Hθ # temporary hamiltonian so master isn't changed
    Hθt = subs(Hθt, the1=>0.0, the2=>0.0, the3=>0.0, the4=>0.0) |> Sym
    Hθtt = Hθt
    Hθtt = subs(Hθtt, A=>x[1,1], B=>x[1,2], S=>x[1,3], W=>x[1,4]) |> float
    xi[1,:] = Δτ*(Hθtt) + x[1,:]
    # End point
    Hθt = subs(Hθt, A=>x[end,1], B=>x[end,2], S=>x[end,3], W=>x[end,4]) |> float
    xi[2,:] = Δτ*(Hθt) + x[end,:]
    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 4)
    for i = 2:NG
        # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
        Hxt = Hx
        Hθθt = Hθθ
        Hθxt = Hθx
        Hxt = subs(Hxt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4]) |> Sym
        Hxt = subs(Hxt, the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
        Hθθt = subs(Hθθt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4]) |> Sym
        Hθθt = subs(Hθθt, the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
        Hθxt = subs(Hθxt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4]) |> Sym
        Hθxt = subs(Hθxt, the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
        for j = 1:4
            K[i,j] = x[i,j]
            for k = 1:4
                K[i,j] -= Δτ*λs[i]*(Hθxt[j,k]*xprim[i,k])
                K[i,j] += Δτ*(Hθθt[j,k]*Hxt[k])
            end
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end

    # Make an initial guess of the path, as prior path
    newxi = x
    # make f! a closure of g! for specific xi, C, K
    f!(F,x) = g!(F,x,C,K,xi,NG)
    # Then put all this into the solver
    newx = nlsolve(f!, newxi)
    return(newx)
end

# main function generates symbolic matrices that the 4 variables can be subbed into
function gMAP(Ω,ϕ,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,Ne,NM,NG,Nmid,Δτ,high2low)
    # make vector of these optimisation parameters
    paras = [ K; k; Q; q; kmin; Kmin; qmin; Qmin; f; r; F ]
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx = gensyms(paras)

    # Now generate an initial path to optimize over
    ss1, sad, ss2 = nullcline(F,r,f,ϕ,K,k,Ne,high2low)
    # a = collect(linspace(ss1[1],ss2[1],NG+1))
    # b = collect(linspace(ss1[2],ss2[2],NG+1))
    # s = collect(linspace(ss1[3],ss2[3],NG+1))
    # w = collect(linspace(ss1[4],ss2[4],NG+1))
    a1 = collect(linspace(ss1[1],sad[1],(NG/2)+1))
    a2 = collect(linspace(sad[1],ss2[1],(NG/2)+1))
    a = vcat(a1,a2[2:length(a2)])
    b1 = collect(linspace(ss1[2],sad[2],(NG/2)+1))
    b2 = collect(linspace(sad[2],ss2[2],(NG/2)+1))
    b = vcat(b1,b2[2:length(b2)])
    s1 = collect(linspace(ss1[3],sad[3],(NG/2)+1))
    s2 = collect(linspace(sad[3],ss2[3],(NG/2)+1))
    s = vcat(s1,s2[2:length(s2)])
    w1 = collect(linspace(ss1[4],sad[4],(NG/2)+1))
    w2 = collect(linspace(sad[4],ss2[4],(NG/2)+1))
    w = vcat(w1,w2[2:length(w2)])
    x = hcat(a,b,s,w)

    # Then appropriatly discretise the path such that it works with this algorithm
    x = discretise(x,NG)
    test = zeros(size(x,1))
    # Set up method to tell if is converged
    convrg = false
    l = 0
    gr()
    while convrg == false
        @time x, xprim, λs, ϑs, λprim = genvars(x,λ,ϑ,NG)
        @time newx = linsys(x,xprim,λs,ϑs,λprim,Hx,Hθ,Hθθ,Hθx,Δτ,NG)
        xn = discretise(newx.zero,NG)
        # delta is the sum of the differences of all the points in the path
        δ = 0
        for i = 1:NG+1
            for j = 1:4
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        S = Ŝ(x,xprim,λs,ϑs,λprim,NG)
        print("$(δ),$(sum(S))\n")
        plot(S)
        savefig("../Results/S$(l).png")
        plot(x[:,3],x[:,4])
        savefig("../Results/SvsW$(l).png")
        l += 1
        # Now overwrite old x
        x = xn
        if δ <=  0.000000005
            convrg = true
            print("$(l) steps to converge\n")
        end
    end
    return(x)
end

# Function to calculate the action of a given path
function Ŝ(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int)
    S = zeros(NG-1)
    for j = 1:4
        S[1] += (3/(2*NG))*(xprim[2,j]*ϑs[2,j])
        # Not excluding the midpoint as the contribution is vanishing
        # Might have to rethink this for the 4 species case
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
            if ts[j] >= t[i] || j == NG+1
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
    path[NM+1,:] = x[NG+1,:]
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

function main()
    # General parameters
    Ω = 300
    ϕ = 0.1 # ratio ϕ = q/k
    K = 10
    k = K*Ω
    Q = K*ϕ
    q = Q*Ω
    kmin = 10.0^-20 # set all too 10.0^-20 for now
    Kmin = 10.0^-20
    qmin = 10.0^-20
    Qmin = 10.0^-20
    f = 1000/(Ω^2) # Promoter switching
    r = 10
    F = 250
    Ne = 12000 # number of elements in the system

    # Optimisation parameters
    NM = 150 # number of segments to discretise MAP onto
    NG = 150 # number of segments to optimize gMAP over
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.0000002#0.0000001 # I've made this choice arbitarily, too large and the algorithm breaks
    high2low = true # Set if starting from high state or low state

    # Now call simulation function with these parameters
    path = gMAP(Ω,ϕ,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,Ne,NM,NG,Nmid,Δτ,high2low)
    x, xprim, λs, ϑs, λprim = genvars(path)
    plot(path[:,1],path[:,2])
    savefig("../Results/Graph1.png")
    plot(path[:,3],path[:,4])
    savefig("../Results/Graph2.png")
end


@time main()
