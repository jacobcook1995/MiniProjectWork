#!/usr/bin/env julia
# gMAP2bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to
# render the path in the usual (and analysable) form

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
using SymEngine
import GR # Need this to stop world age plotting error

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function nullcline(ps::Array{Float64,1},high2low::Bool)
    # gonna use SymPy here
    # A1(x) = sqrt.((r/f)*(q./((qmin+Q).*x-Qmin) - 1))
    # A2(x) = (1/(kmin+K))*((k*r)./(r+f*x.^2) + Kmin)
    A1(x) = real(sqrt(complex((ps[10]/ps[9])*(ps[4]/((ps[3]+ps[7])*x - ps[8]) - 1))))
    A2(x) = (1/(ps[5]+ps[1]))*((ps[2]*ps[10])/(ps[10]+ps[9]*x^2) + ps[6])
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    while three == false
        bs = fzeros(g, 0, 15.0)
        n = length(bs)
        gs = 0
        bad = zeros(Int64,0)
        for i = 1:n
            # check if this is actual solution and not artifact of the numerical method
            gs = g(bs[i])
            tol = 1.0e-14
            if gs >= 0 + tol || gs <= 0 - tol
                bad = append!(bad,i)
            end
        end
        if length(bad) != 0
            n = n - length(bad)
            bs = deleteat!(bs, bad)
        end
        if n == 3
            three = true
        end
    end
    sad = [ A1(bs[2]), bs[2] ]
    if high2low == true
        ss1 = [ A1(bs[1]), bs[1] ]
        ss2 = [ A1(bs[3]), bs[3] ]
    else
        ss1 = [ A1(bs[3]), bs[3] ]
        ss2 = [ A1(bs[1]), bs[1] ]
    end
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

# make a symbolic diffusion matrix
function Ds()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,2,2)
    e[1,2] = e[2,1] = 0
    e[1,1] = sqrt(k*r/(r + f*B^2) + kmin*A + K*A + Kmin)
    e[2,2] = sqrt(q*r/(r + f*A^2) + qmin*B + Q*B + Qmin)
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    Dmin = inv(D)
    for j = 1:2
        for i = 1:2
            Dmin[i,j] = expand(Dmin[i,j])
        end
    end
    return(Dmin)
end

# function to make a symbolic equation vector
function bs()
    A, B, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("A B K k Q q kmin Kmin qmin Qmin f r")
    # Make a symbolic version of the matrix, needs no input in this case
    b = Array{SymEngine.Basic,1}(undef,2)
    b[1] = k*r/(r + f*B^2) - kmin*A - K*A + Kmin
    b[2] = q*r/(r + f*A^2) - qmin*B - Q*B + Qmin
    return(b)
end

# function to generate a symbolic equation for the Hamiltonian at point (X, θ)
function Hs()
    the1, the2 = symbols("the1 the2")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the1*b[1] + the2*b[2]
    H += 0.5*the1*(D[1,1]*the1 + D[1,2]*the2)
    H += 0.5*the2*(D[2,1]*the1 + D[2,2]*the2)
    H = expand(H)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    A, B = symbols("A B")
    # generate Hamiltonian
    H = Hs()
    Hx = Array{SymEngine.Basic,1}(undef,2)
    Hx[1] = diff(H, A)
    Hx[2] = diff(H, B)
    return(Hx)
end

# function to generate first differential of the symbolic hamiltonian in θ
function Hθs()
    the1, the2 = symbols("the1 the2")
    # generate Hamiltonian
    H = Hs()
    Hθ = Array{SymEngine.Basic,1}(undef,2)
    Hθ[1] = diff(H, the1)
    Hθ[2] = diff(H, the2)
    return(Hθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ
function Hθθs()
    the1, the2 = symbols("the1 the2")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθθ = Array{SymEngine.Basic,2}(undef,2,2)
    Hθθ[1,1] = diff(Hθ[1],the1)
    Hθθ[1,2] = diff(Hθ[1],the2)
    Hθθ[2,1] = diff(Hθ[2],the1)
    Hθθ[2,2] = diff(Hθ[2],the2)
    return(Hθθ)
end

# function to generate the second differential of the symbolic hamiltonian in θ follwed by X
function Hθxs()
    A, B = symbols("A B")
    # generate Hamiltonian
    Hθ = Hθs()
    Hθx = Array{SymEngine.Basic,2}(undef,2,2)
    Hθx[1,1] = diff(Hθ[1],A)
    Hθx[1,2] = diff(Hθ[1],B)
    Hθx[2,1] = diff(Hθ[2],A)
    Hθx[2,2] = diff(Hθ[2],B)
    return(Hθx)
end

# function to find a symbolic equation for λ the determenistic speed
function λs()
    y1, y2 = symbols("y1 y2")
    b = bs()
    Dmin = Dmins()
    num = 0
    num += b[1]*Dmin[1,1]*b[1] + b[1]*Dmin[1,2]*b[2]
    num += b[2]*Dmin[2,1]*b[1] + b[2]*Dmin[2,2]*b[2]
    den = 0
    den += y1*Dmin[1,1]*y1 + y1*Dmin[1,2]*y2
    den += y2*Dmin[2,1]*y1 + y2*Dmin[2,2]*y2
    λ = sqrt(num)/sqrt(den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2 = symbols("y1 y2")
    λ = λs()
    Dmin = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,2)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    ϑ = Array{SymEngine.Basic,1}(undef,2)
    ϑ[1] = Dmin[1,1]*c[1] + Dmin[1,2]*c[2]
    ϑ[2] = Dmin[2,1]*c[1] + Dmin[2,2]*c[2]
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
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = symbols("K k Q q kmin Kmin qmin Qmin f r")
    # now perform substitutions
    # ps not right here CHANGE! CHNAGE! CHANGE!
    λ = subs(λ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
    λ = subs(λ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
    H = subs(H, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
    H = subs(H, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
    for i = 1:2
        ϑ[i] = subs(ϑ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        ϑ[i] = subs(ϑ[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
        Hθ[i] = subs(Hθ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        Hθ[i] = subs(Hθ[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
        Hx[i] = subs(Hx[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        Hx[i] = subs(Hx[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
    end
    for j = 1:2
        for i = 1:2
            Hθx[i,j] = subs(Hθx[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
            Hθx[i,j] = subs(Hθx[i,j], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
            Hθθ[i,j] = subs(Hθθ[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
            Hθθ[i,j] = subs(Hθθ[i,j], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10])
        end
    end
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to generate the variables needed for a given algoritm iteration
function genvars(x::AbstractArray,λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int,Nmid::Int)
    # define neccesary symbols
    A, B, y1, y2 = symbols("A B y1 y2")
    # calculate velocities
    xprim = fill(NaN, NG+1, 2)
    for i = 2:NG
        for j = 1:2
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:Nmid-1
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = subs(λt, A=>x[i,1], B=>x[i,2])
        λs[i] = subs(λt, y1=>xprim[i,1], y2=>xprim[i,2]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 2)
    ϑt = Array{SymEngine.Basic,1}(undef,2)
    for j = 1:2
        for i = 2:NG
            ϑt[j] = subs(ϑ[j], A=>x[i,1], B=>x[i,2])
            ϑs[i,j] = subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2]) |> float
        end
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end
    return(x,xprim,λs,ϑs,λprim)
end

# function to be solved by NLsolve
function g!(F::AbstractArray,x::AbstractArray,C::AbstractVector,K::AbstractArray,xi::AbstractArray,NG::Int64,Nmid::Int64)
    # Start point
    F[1,1] = x[1,1] - xi[1,1]
    F[1,2] = x[1,2] - xi[1,2]
    # first path
    for i = 2:Nmid-1
        for j = 1:2
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # midpoint
    F[Nmid,1] = x[Nmid,1] - xi[2,1]
    F[Nmid,2] = x[Nmid,2] - xi[2,2]
    # second path
    for i = Nmid+1:NG
        for j = 1:2
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # end point
    F[NG+1,1] = x[NG+1,1] - xi[3,1]
    F[NG+1,2] = x[NG+1,2] - xi[3,2]
    return(F)
end

# function to solve the system of linear equations
function linsys(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,
                Hx::Array{SymEngine.Basic,1},Hθ::Array{SymEngine.Basic,1},Hθθ::Array{SymEngine.Basic,2},
                Hθx::Array{SymEngine.Basic,2},Δτ::Number,NG::Int,Nmid::Int,H::SymEngine.Basic)
    # define relevant symbols
    A, B, the1, the2 = symbols("A B the1 the2")
    # Make array to store fixed points
    xi = fill(NaN, 3, 2)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Array{SymEngine.Basic,1}(undef,2) # temporary hamiltonian so master isn't changed
    Hxt = Array{SymEngine.Basic,1}(undef,2)
    for i = 1:2
        Hθt[i] = subs(Hθ[i], the1=>0.0, the2=>0.0)
        Hxt[i] = subs(Hx[i], the1=>0.0, the2=>0.0)
    end
    Hθtt = Array{SymEngine.Basic,1}(undef,2)
    Hθttt = Array{SymEngine.Basic,1}(undef,2)
    Ht = subs(H, the1=>0.0, the2=>0.0)
    Ht = subs(Ht, A=>x[Nmid,1], B=>x[Nmid,2]) |> float
    # loop to define fixed points
    for i = 1:2
        Hxt[i] = subs(Hxt[i], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        Hθttt[i] = subs(Hθt[i], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        Hθtt[i] = subs(Hθt[i], A=>x[1,1], B=>x[1,2]) |> float
        Hθt[i] = subs(Hθt[i], A=>x[end,1], B=>x[end,2]) |> float
        # start point
        xi[1,i] = Δτ*(Hθtt[i]) + x[1,i]
        # midpoint
        xi[2,i] = -Δτ*(Hxt[i]*Ht) + x[Nmid,i]
        # End point
        xi[3,i] = Δτ*(Hθt[i]) + x[end,i]
    end
    Hxθt = Array{SymEngine.Basic,2}(undef,2,2)
    # loop to calculate Hxθ
    for j = 1:2
        for i = 1:2
            Hxθt[i,j] = subs(Hθx[i,j], the1=>0.0, the2=>0.0)
            Hxθt[i,j] = subs(Hxθt[i,j], A=>x[Nmid,1], B=>x[Nmid,2]) |> float
        end
    end
    # now transpose to get right form
    Hxθt = transpose(Hxθt) # possible source of error if ive misunderstood here
    # loop for additional midpoint terms
    for j = 1:2
        for i = 1:2
            xi[2,i] -= Δτ*(Hxθt[i,j]*Hθt[j])
        end
    end

    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 2)
    Hxt = Array{SymEngine.Basic,1}(undef,2)
    Hθθt = Array{SymEngine.Basic,2}(undef,2,2)
    Hθxt = Array{SymEngine.Basic,2}(undef,2,2)
    # Put initial values for K in
    for j = 1:2
        for i = 2:NG
            K[i,j] = x[i,j]
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end
    for l = 1:2
        for i = 2:NG
            # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
            Hxt[l] = subs(Hx[l], A=>x[i,1], B=>x[i,2])
            Hxt[l] = subs(Hxt[l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
            for m = 1:2
                Hθθt[m,l] = subs(Hθθ[m,l], A=>x[i,1], B=>x[i,2])
                Hθθt[m,l] = subs(Hθθt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
                Hθxt[m,l] = subs(Hθx[m,l], A=>x[i,1], B=>x[i,2])
                Hθxt[m,l] = subs(Hθxt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2]) |> float
                # Update K's with new contributions from Hamiltonians
                K[i,m] -= Δτ*λs[i]*(Hθxt[m,l]*xprim[i,l])
                K[i,m] += Δτ*(Hθθt[m,l]*Hxt[l])
            end
        end
    end
    # Make an initial guess of the path, as prior path
    newxi = x
    # make f! a closure of g! for specific xi, C, K
    f!(F,x) = g!(F,x,C,K,xi,NG,Nmid)
    # Then put all this into the solver
    newx = nlsolve(f!, newxi)
    return(newx)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray,NG::Int64,Nmid::Int64,spline::Bool,sub::Int64)
    if spline == false
        # need to interpolate the data from x onto disx, preserving |x'| = const,
        # i.e equal distance between points
        s1 = zeros(Nmid)
        s2 = zeros(NG+2-Nmid)
        s1[1] = 0
        s2[1] = 0
        for i = 2:Nmid
            dA = x[i,1] - x[i-1,1]
            dB = x[i,2] - x[i-1,2]
            s1[i] = s1[i-1] + sqrt(dA^2 + dB^2) # Could probably drop the sqrts to speed up the code
        end
        for i = 2:(NG+2-Nmid)
            dA = x[i+Nmid-1,1] - x[i+Nmid-2,1]
            dB = x[i+Nmid-1,2] - x[i+Nmid-2,2]
            s2[i] = s2[i-1] + sqrt(dA^2 + dB^2) # Could probably drop the sqrts to speed up the code
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
        disx = zeros(NG+1,2)
        disx[1,:] = x[1,:]
        disx[Nmid,:] = x[Nmid,:]
        disx[NG+1,:] = x[NG+1,:]
        # This is done to linear order, which is probably good enough
        for i = 2:Nmid-1
            one = inds1[i] - 1
            two = inds1[i]
            s₀ = s1[one]
            s₁ = s1[two]
            for j = 1:2
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
            for j = 1:2
                x₀ = x[one,j]
                x₁ = x[two,j]
                disx[i,j] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
            end
        end
    else
        # make spline 1
        A1 = zeros(Nmid,Nmid)
        b1 = zeros(Nmid)
        for i = 1:Nmid
            if i == 1
                b1[i] = 3*(x[i+1,2] - x[i,2])/((x[i+1,1] - x[i,1])^2)
                A1[i,i] = 2/(x[i+1,1] - x[i,1])
                A1[i,i+1] = 1/(x[i+1,1] - x[i,1])
            elseif i == Nmid
                b1[i] = 3*(x[i,2] - x[i-1,2])/((x[i,1] - x[i-1,1])^2)
                A1[i,i-1] = 1/(x[i,1] - x[i-1,1])
                A1[i,i] = 2/(x[i,1] - x[i-1,1])
            else
                b1[i] = 3*(x[i+1,2] - x[i,2])/((x[i+1,1] - x[i,1])^2) + 3*(x[i,2] - x[i-1,2])/((x[i,1] - x[i-1,1])^2)
                A1[i,i-1] = 1/(x[i,1] - x[i-1,1])
                A1[i,i] = 2*(1/(x[i,1] - x[i-1,1]) + 1/(x[i+1,1] - x[i,1]))
                A1[i,i+1] = 1/(x[i+1,1] - x[i,1])
            end
        end
        # make spline 2
        A2 = zeros(NG+2-Nmid,NG+2-Nmid)
        b2 = zeros(NG+2-Nmid)
        for i = Nmid:(NG+1)
            j = i - Nmid + 1
            if i == Nmid
                b2[j] = 3*(x[i+1,2] - x[i,2])/((x[i+1,1] - x[i,1])^2)
                A2[j,j] = 2/(x[i+1,1] - x[i,1])
                A2[j,j+1] = 1/(x[i+1,1] - x[i,1])
            elseif i == NG+1
                b2[j] = 3*(x[i,2] - x[i-1,2])/((x[i,1] - x[i-1,1])^2)
                A2[j,j-1] = 1/(x[i,1] - x[i-1,1])
                A2[j,j] = 2/(x[i,1] - x[i-1,1])
            else
                b2[j] = 3*(x[i+1,2] - x[i,2])/((x[i+1,1] - x[i,1])^2) + 3*(x[i,2] - x[i-1,2])/((x[i,1] - x[i-1,1])^2)
                A2[j,j-1] = 1/(x[i,1] - x[i-1,1])
                A2[j,j] = 2*(1/(x[i,1] - x[i-1,1]) + 1/(x[i+1,1] - x[i,1]))
                A2[j,j+1] = 1/(x[i+1,1] - x[i,1])
            end
        end
        # now calculate inverse matrices
        A1 = inv(A1)
        A2 = inv(A2)
        k1 = A1*b1
        k2 = A2*b2
        # need to interpolate the data from x onto disx, preserving |x'| = const,
        # i.e equal distance between points
        n1 = Nmid + (sub-1)*(Nmid-1)
        s1 = zeros(n1)
        x1 = zeros(n1)
        y1 = zeros(n1)
        s1[1] = 0
        x1[1] = x[1,1]
        y1[1] = x[1,2]
        j = 1
        # determine a and b values for initial segement
        a = k1[j]*(x[j+1,1] - x[j,1]) - (x[j+1,2] - x[j,2])
        b = -k1[j+1]*(x[j+1,1] - x[j,1]) + (x[j+1,2] - x[j,2])
        for i = 2:n1
            x1[i] = x1[i-1] + (1/sub)*(x[j+1,1] - x[j,1])
            t = (x1[i] - x[j,1])/(x[j+1,1] - x[j,1])
            y1[i] = (1 - t)*x[j,2] + t*x[j+1,2] + t*(1-t)*(a*(1-t) + b*t)
            s1[i] = s1[i-1] + sqrt((x1[i]-x1[i-1])^2 + (y1[i]-y1[i-1])^2)
            if (i-1) % sub == 0 && i != n1
                j += 1
                # redefine when moving onto another segment
                a = k1[j]*(x[j+1,1] - x[j,1]) - (x[j+1,2] - x[j,2])
                b = -k1[j+1]*(x[j+1,1] - x[j,1]) + (x[j+1,2] - x[j,2])
            end
        end
        n2 = NG + 2 - Nmid + (sub-1)*(NG+1-Nmid)
        s2 = zeros(n2)
        x2 = zeros(n2)
        y2 = zeros(n2)
        s2[1] = 0
        x2[1] = x[Nmid,1]
        y2[1] = x[Nmid,2]
        j = Nmid
        # determine a and b values for initial segement
        a = k2[j+1-Nmid]*(x[j+1,1] - x[j,1]) - (x[j+1,2] - x[j,2])
        b = -k2[j+2-Nmid]*(x[j+1,1] - x[j,1]) + (x[j+1,2] - x[j,2])
        for i = 2:n2
            x2[i] = x2[i-1] + (1/sub)*(x[j+1,1] - x[j,1])
            t = (x2[i] - x[j,1])/(x[j+1,1] - x[j,1])
            y2[i] = (1 - t)*x[j,2] + t*x[j+1,2] + t*(1-t)*(a*(1-t) + b*t)
            s2[i] = s2[i-1] + sqrt((x2[i]-x2[i-1])^2 + (y2[i]-y2[i-1])^2)
            if (i-1) % sub == 0 && i != n2
                j += 1
                # redefine when moving onto another segment
                a = k2[j+1-Nmid]*(x[j+1,1] - x[j,1]) - (x[j+1,2] - x[j,2])
                b = -k2[j+2-Nmid]*(x[j+1,1] - x[j,1]) + (x[j+1,2] - x[j,2])
            end
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
                if s1[j] >= ls1[i] || j == n1
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
                if s2[j] >= ls2[i] || j == n2
                    inds2[i] = j + n1 - 1
                    higher = true
                else
                    j += 1
                end
            end
        end
        # First do mid points and end points as they should be fixed
        disx = zeros(NG+1,2)
        disx[1,:] = x[1,:]
        disx[Nmid,:] = x[Nmid,:]
        disx[NG+1,:] = x[NG+1,:]
        # This is done to linear order, which is probably good enough
        for i = 2:Nmid-1
            one = inds1[i] - 1
            two = inds1[i]
            s₀ = s1[one]
            s₁ = s1[two]
            x₀ = x1[one]
            x₁ = x1[two]
            disx[i,1] = x₀ + (ls1[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
            y₀ = y1[one]
            y₁ = y1[two]
            disx[i,2] = y₀ + (ls1[i] - s₀)*(y₁ - y₀)/(s₁ - s₀)
        end

        for i = Nmid+1:NG
            one = inds2[i+1-Nmid] - 1
            two = inds2[i+1-Nmid]
            s₀ = s2[one+1-n1]
            s₁ = s2[two+1-n1]
            x₀ = x2[one+1-n1]
            x₁ = x2[two+1-n1]
            disx[i,1] = x₀ + (ls2[i+1-Nmid] - s₀)*(x₁ - x₀)/(s₁ - s₀)
            y₀ = y2[one+1-n1]
            y₁ = y2[two+1-n1]
            disx[i,2] = y₀ + (ls2[i+1-Nmid] - s₀)*(y₁ - y₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

# Function to generate an initial path then run the alorithm until a MAP is obtained
# This path is then returned for other functions
function gMAP(ps::Array{Float64,1},NG::Int64,Nmid::Int64,Δτ::Float64,high2low::Bool,spline::Bool)
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(ps)
    # First find the steady states and saddle point
    ss1, sad, ss2 = nullcline(ps,high2low)
    len = round(Int64,(NG/2)+1) # must be integer
    a1 = collect(range(ss1[1],stop=sad[1],length=len))
    a2 = collect(range(sad[1],stop=ss2[1],length=len))
    a = vcat(a1,a2[2:end])
    b1 = collect(range(ss1[2],stop=sad[2],length=len))
    b2 = collect(range(sad[2],stop=ss2[2],length=len))
    b = vcat(b1,b2[2:end])
    x = hcat(a,b)
    # Then appropriatly discretise the path such that it works with this algorithm
    prec = 1000
    x = discretise(x,NG,Nmid,spline,prec)
    # Set up method to tell if is converged
    convrg = false
    l = 0
    while convrg == false
        x, xprim, λs, ϑs, λprim = genvars(x,λ,ϑ,NG,Nmid)
        newx = linsys(x,xprim,λs,ϑs,λprim,Hx,Hθ,Hθθ,Hθx,Δτ,NG,Nmid,H)
        xn = discretise(newx.zero,NG,Nmid,spline,prec)
        # delta is the sum of the differences of all the points in the path
        δ = 0
        for i = 1:NG+1
            for j = 1:2
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        l += 1
        if l % 200 == 0
            println("$l,$δ")
        end
        # Now overwrite old x
        x = xn
        if δ <= 0.000000005#0.00000000005
            convrg = true
            print("$(l) steps to converge\n")
        end
    end
    return(x,λ,ϑ)
end

# Function to calculate the action of a given path
function Ŝ(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int64)
    S = zeros(NG-1)
    S[1] = (3/(2*NG))*(xprim[2,1]*ϑs[2,1] + xprim[2,2]*ϑs[2,2])
    # Not excluding the midpoint as the contribution is vanishing
    # Might have to rethink this for the 4 species case
    for i = 3:NG-1
        S[i-1] += (1/NG)*(xprim[i,1]*ϑs[i,1] + xprim[i,2]*ϑs[i,2])
    end
    S[NG-1] += (3/(2*NG))*(xprim[NG,1]*ϑs[NG,1] + xprim[NG,2]*ϑs[NG,2])
    return(S)
end

# function to find the times of each point
function times(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,NG::Int64)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::AbstractVector,x::AbstractArray,NG::Int64,NM::Int64)
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
    path = zeros(NM+1,2)
    path[1,:] = x[1,:]
    path[NM+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        for j = 1:2
            x₀ = x[one,j]
            x₁ = x[two,j]
            path[i,j] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
        end
    end
    return(path)
end

function params(Ah::BigFloat,Al::BigFloat,Bh::BigFloat,Bl::BigFloat,F::BigFloat)
    Q = Qmin = K = Kmin = k = kmin = qmin = q = r = f = convert(BigFloat,0.0)
    Q = F/(Bh - Bl)
    Qmin = Q*Bl
    K = F/(Ah - Al)
    Kmin = K*Al
    correct = false
    l = 0
    while correct == false
        r = 10000*rand()
        f = 10000*rand()
        r1 = ((1/(r + f*Bl^2))-(1/(r + f*Bh^2)))/((1/(r + f*Bl^2))+(1/(r + f*Bh^2)))
        kmin = F*(1-r1)/((Ah+Al)*r1 - Ah + Al)
        k = (F + kmin*(Ah+Al))/(r*(1/(r + f*Bl^2) + 1/(r + f*Bh^2)))
        r2 = ((1/(r + f*Al^2))-(1/(r + f*Ah^2)))/((1/(r + f*Al^2))+(1/(r + f*Ah^2)))
        qmin = F*(1-r2)/((Bh+Bl)*r2 - Bh + Bl)
        q = (F + qmin*(Bh+Bl))/(r*(1/(r + f*Al^2) + 1/(r + f*Ah^2)))
        if kmin > 0 && k > 0 && qmin > 0 && q > 0
            correct = true
            println("$(l+1) attempts to generate parameters")
        else
            l += 1
        end
    end
    return(K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r)
end

function main()
    # Firstly should generate parameters
    # define flux over high state
    F = convert(BigFloat,10.0)
    Ah = convert(BigFloat,10.0)
    Al = convert(BigFloat,1.0)
    Bh = convert(BigFloat,10.1)
    Bl = convert(BigFloat,2.0)

    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r = params(Ah,Al,Bh,Bl,F)
    # put in vector
    ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
    ps = convert(Array{Float64,1},ps)

    # Set parameters of the optimization
    NM = 600 # number of segments to discretise MAP onto
    NG = 600 # number of segments to optimize gMAP over
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.005#0.001 # I've made this choice arbitarily, too large and the algorithm breaks
    high2low = false
    spline = false

    path, λ, ϑ = gMAP(ps,NG,Nmid,Δτ,high2low,spline)
    x, xprim, λs, ϑs, λprim = genvars(path,λ,ϑ,NG,Nmid)
    # use function Ŝ to find the action associated with this path
    λs[1] = λs[2]
    λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
    λs[end] = λs[end-1]
    S = Ŝ(x,xprim,λs,ϑs,λprim,NG)
    print("Associated Action = $(sum(S))\n")
    tims1 = times(x,xprim,λs,ϑs,λprim,NG)
    print("Time of path = $(tims1[end])\n")
    path1 = timdis(tims1,x,NG,NM)
    # now run reverse path
    path, λ, ϑ = gMAP(ps,NG,Nmid,Δτ,~high2low,spline)
    x, xprim, λs, ϑs, λprim = genvars(path,λ,ϑ,NG,Nmid)
    # use function Ŝ to find the action associated with this path
    λs[1] = λs[2]
    λs[Nmid] = (λs[Nmid+1] + λs[Nmid-1])/2
    λs[end] = λs[end-1]
    S = Ŝ(x,xprim,λs,ϑs,λprim,NG)
    print("Associated Action = $(sum(S))\n")
    tims2 = times(x,xprim,λs,ϑs,λprim,NG)
    print("Time of path = $(tims2[end])\n")
    path2 = timdis(tims2,x,NG,NM)

    # Block of code to write all this data to a file so I can go through it
    if length(ARGS) >= 1
        output_file1 = "../Results/0708/$(ARGS[1])1.csv"
        out_file1 = open(output_file1, "w")
        # open file for writing
        for i = 1:size(path1,1)
            line = "$(path1[i,1]),$(path1[i,2])\n"
            write(out_file1, line)
        end
        close(out_file1)
        output_file2 = "../Results/0708/$(ARGS[1])2.csv"
        out_file2 = open(output_file2, "w")
        # open file for writing
        for i = 1:size(path2,1)
            line = "$(path2[i,1]),$(path2[i,2])\n"
            write(out_file2, line)
        end
        close(out_file2)
        output_filep = "../Results/0708/$(ARGS[1])p.csv"
        out_filep = open(output_filep, "w")
        # open file for writing
        for i = 1:size(ps,1)
            line = "$(ps[i])\n"
            write(out_filep, line)
        end
        line = "$(tims1[end])\n"
        write(out_filep, line)
        line = "$(tims2[end])\n"
        write(out_filep, line)
        close(out_filep)
    end
end

@time main()
