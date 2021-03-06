#!/usr/bin/env julia
# gMAP4bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to this
# time done for the 4 species case

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
using SymEngine # Trying alternative substitution method to see if this is faster
using SymPy
using LinearAlgebra
ENV["GKSwstype"]="png" # ????????? this might stop the GKSterm error
import GR # Need this to stop world age plotting error?

# first should make generic functions that take in a full set of parameters, that
# can perform the machine alegbra and then sub in the parameter set to generate
# numeric values for each of the instances considered

# make a symbolic diffusion matrix
function Ds()
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = SymEngine.symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
    # Make a symbolic version of the matrix, needs no input in this case
    e = Array{SymEngine.Basic,2}(undef,4,9)
    e[1,1] = sqrt(k*S*r/(r + f*B^2))
    e[1,2] = -sqrt(kmin*A)
    e[1,3] = 0
    e[1,4] = 0
    e[1,5] = -sqrt(K*A)
    e[1,6] = sqrt(Kmin*W)
    e[1,7] = 0
    e[1,8] = 0
    e[1,9] = 0
    e[2,1] = 0
    e[2,2] = 0
    e[2,3] = sqrt(q*S*r/(r + f*A^2))
    e[2,4] = -sqrt(qmin*B)
    e[2,5] = 0
    e[2,6] = 0
    e[2,7] = -sqrt(Q*B)
    e[2,8] = sqrt(Qmin*W)
    e[2,9] = 0
    e[3,1] = -sqrt(k*S*r/(r + f*B^2))
    e[3,2] = sqrt(kmin*A)
    e[3,3] = -sqrt(q*S*r/(r + f*A^2))
    e[3,4] = sqrt(qmin*B)
    e[3,5] = 0
    e[3,6] = 0
    e[3,7] = 0
    e[3,8] = 0
    e[3,9] = 0.001# peak pretty unmanagable at (0.0001) really small spike at (0.01) but totals start to change substaintailly
    e[4,1] = 0
    e[4,2] = 0
    e[4,3] = 0
    e[4,4] = 0
    e[4,5] = sqrt(K*A)
    e[4,6] = -sqrt(Kmin*W)
    e[4,7] = sqrt(Q*B)
    e[4,8] = -sqrt(Qmin*W)
    e[4,9] = 0
    # Now do the transformations required
    eT = transpose(e)
    D = e*eT
    return(D)
end

# function to generate symbolic inverse diffusion matrix
function Dmins()
    D = Ds()
    detD = det4x4(D)
    cofacD = cofac4x4(D)
    adjD = transpose(cofacD)
    return(detD,adjD)
end

# function to return determinent of a 2X2 matrix
function det2x2(Mat::Array{SymEngine.Basic,2})
    if size(Mat,1) != 2 || size(Mat,2) != 2
        println("Error:Wrong matrix size must be 2x2")
        error()
    else
        Det = Mat[1,1]*Mat[2,2] - Mat[1,2]*Mat[2,1]
        return(Det)
    end
end

# function to return determinent of a 3X3 matrix
function det3x3(Mat::Array{SymEngine.Basic,2})
    if size(Mat,1) != 3 || size(Mat,2) != 3
        println("Error:Wrong matrix size must be 3x3")
        error()
    else
        Det = Mat[1,1]*det2x2(Mat[2:3,2:3]) - Mat[1,2]*det2x2(Mat[2:3,1:2:3]) + Mat[1,3]*det2x2(Mat[2:3,1:2])
        return(Det)
    end
end

# function to return determinent of a 4X4 matrix
function det4x4(Mat::Array{SymEngine.Basic,2})
    if size(Mat,1) != 4 || size(Mat,2) != 4
        println("Error:Wrong matrix size must be 4x4")
        error()
    else
        Det = Mat[1,1]*det3x3(Mat[2:4,2:4]) - Mat[1,2]*det3x3(Mat[2:4,[1,3,4]]) + Mat[1,3]*det3x3(Mat[2:4,[1,2,4]]) - Mat[1,4]*det3x3(Mat[2:4,1:3])
        return(Det)
    end
end

# function to return cofactor matrix for a 4X4 matrix
function cofac4x4(Mat::Array{SymEngine.Basic,2})
    if size(Mat,1) != 4 || size(Mat,2) != 4
        println("Error:Wrong matrix size must be 4x4")
        error()
    else
        Cofac = Array{SymEngine.Basic,2}(undef,4,4)
        for i = 1:4
            ilist = setdiff([1,2,3,4], [i])
            for j = 1:4
                jlist = setdiff([1,2,3,4], [j])
                Cofac[i,j] = ((-1)^(i+j))*det3x3(Mat[ilist,jlist])
            end
        end
        return(Cofac)
    end
end

# function to make a symbolic equation vector
function bs()
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = SymEngine.symbols("A B S W K k Q q kmin Kmin qmin Qmin f r F")
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
    the1, the2, the3, the4 = SymEngine.symbols("the1 the2 the3 the4")
    # generate symbolic arrays for b and D
    b = bs()
    D = Ds()
    H = 0
    H += the1*b[1] + the2*b[2] + the3*b[3] + the4*b[4]
    H += 0.5*the1*(D[1,1]*the1 + D[1,2]*the2 + D[1,3]*the3 + D[1,4]*the4)
    H += 0.5*the2*(D[2,1]*the1 + D[2,2]*the2 + D[2,3]*the3 + D[2,4]*the4)
    H += 0.5*the3*(D[3,1]*the1 + D[3,2]*the2 + D[3,3]*the3 + D[3,4]*the4)
    H += 0.5*the4*(D[4,1]*the1 + D[4,2]*the2 + D[4,3]*the3 + D[4,4]*the4)
    H = SymEngine.expand(H)
    return(H)
end

# function to generate first differential of the symbolic hamiltonian in x
function Hxs()
    A, B, S, W = SymEngine.symbols("A B S W")
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
    the1, the2, the3, the4 = SymEngine.symbols("the1 the2 the3 the4")
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
    the1, the2, the3, the4 = SymEngine.symbols("the1 the2 the3 the4")
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
    A, B, S, W = SymEngine.symbols("A B S W")
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
    y1, y2, y3, y4 = SymEngine.symbols("y1 y2 y3 y4")
    b = bs()
    _, adjD = Dmins()
    num = 0
    num += b[1]*adjD[1,1]*b[1] + b[1]*adjD[1,2]*b[2] + b[1]*adjD[1,3]*b[3] + b[1]*adjD[1,4]*b[4]
    num += b[2]*adjD[2,1]*b[1] + b[2]*adjD[2,2]*b[2] + b[2]*adjD[2,3]*b[3] + b[2]*adjD[2,4]*b[4]
    num += b[3]*adjD[3,1]*b[1] + b[3]*adjD[3,2]*b[2] + b[3]*adjD[3,3]*b[3] + b[3]*adjD[3,4]*b[4]
    num += b[4]*adjD[4,1]*b[1] + b[4]*adjD[4,2]*b[2] + b[4]*adjD[4,3]*b[3] + b[4]*adjD[4,4]*b[4]
    den = 0
    den += y1*adjD[1,1]*y1 + y1*adjD[1,2]*y2 + y1*adjD[1,3]*y3 + y1*adjD[1,4]*y4
    den += y2*adjD[2,1]*y1 + y2*adjD[2,2]*y2 + y2*adjD[2,3]*y3 + y2*adjD[2,4]*y4
    den += y3*adjD[3,1]*y1 + y3*adjD[3,2]*y2 + y3*adjD[3,3]*y3 + y3*adjD[3,4]*y4
    den += y4*adjD[4,1]*y1 + y4*adjD[4,2]*y2 + y4*adjD[4,3]*y3 + y4*adjD[4,4]*y4
    λ = sqrt(num/den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2, y3, y4 = SymEngine.symbols("y1 y2 y3 y4")
    λ = λs()
    detD, adjD = Dmins()
    b = bs()
    c = Array{SymEngine.Basic,1}(undef,4)
    c[1] = λ*y1 - b[1]
    c[2] = λ*y2 - b[2]
    c[3] = λ*y3 - b[3]
    c[4] = λ*y4 - b[4]
    ϑ = Array{SymEngine.Basic,1}(undef,4)
    # theres an additional simplification that could be acheived by multiplying all ϑ by detD
    ϑ[1] = (adjD[1,1]*c[1] + adjD[1,2]*c[2] + adjD[1,3]*c[3] + adjD[1,4]*c[4])
    ϑ[2] = (adjD[2,1]*c[1] + adjD[2,2]*c[2] + adjD[2,3]*c[3] + adjD[2,4]*c[4])
    ϑ[3] = (adjD[3,1]*c[1] + adjD[3,2]*c[2] + adjD[3,3]*c[3] + adjD[3,4]*c[4])
    ϑ[4] = (adjD[4,1]*c[1] + adjD[4,2]*c[2] + adjD[4,3]*c[3] + adjD[4,4]*c[4])
    ϑ = ϑ/detD
    return(ϑ)
end

# function to generate one of each symbolic object and to substitute parameter values into it
# this is a computationally intensive function, need to minimize the number of calls
function gensyms(ps::Array{Float64,1})
    # create symbolic objects
    ϑ = ϑs()
    λ = λs()
    Hθ = Hθs()
    Hθx = Hθxs()
    Hθθ = Hθθs()
    Hx = Hxs()
    H = Hs()
    # specify symbols that will be substituted for
    K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F = SymEngine.symbols("K k Q q kmin Kmin qmin Qmin f r F")
    # now perform substitutions
    λ = SymEngine.subs(λ, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
    λ = SymEngine.subs(λ, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
    H = SymEngine.subs(H, K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
    H = SymEngine.subs(H, qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
    for i = 1:4
        ϑ[i] = SymEngine.subs(ϑ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        ϑ[i] = SymEngine.subs(ϑ[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
        Hθ[i] = SymEngine.subs(Hθ[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        Hθ[i] = SymEngine.subs(Hθ[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
        Hx[i] = SymEngine.subs(Hx[i], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
        Hx[i] = SymEngine.subs(Hx[i], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
    end
    for j = 1:4
        for i = 1:4
            Hθx[i,j] = SymEngine.subs(Hθx[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
            Hθx[i,j] = SymEngine.subs(Hθx[i,j], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
            Hθθ[i,j] = SymEngine.subs(Hθθ[i,j], K=>ps[1], k=>ps[2], Q=>ps[3], q=>ps[4], kmin=>ps[5], Kmin=>ps[6])
            Hθθ[i,j] = SymEngine.subs(Hθθ[i,j], qmin=>ps[7], Qmin=>ps[8], f=>ps[9], r=>ps[10], F=>ps[11])
        end
    end
    return(ϑ,λ,Hθ,Hθx,Hθθ,Hx,H)
end

# function to find nullclines
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N ]
# maxr = range to look over, δ step size in this range
function nullcline(ps::Array{Float64,1},maxr::Float64,δ::Float64)
    # define symbols
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, N = SymPy.symbols("A,B,S,W,K k,Q,q,kmin,Kmin,qmin,Qmin,f,r,F,N")
    # equation for S
    Se = N - W - A - B
    # equation for W
    We = (K*A + Q*B - F)/(Qmin + Kmin)
    # eliminate W from equation for S
    Se = SymPy.subs(Se,W=>We)
    # equation for dB/dt
    dB = (q*S*r)/(r + f*A^2) - qmin*B - Q*B + Qmin*W
    # sub in S and W expressions
    dB = SymPy.subs(dB,S=>Se,W=>We)
    # solve this expression for B
    Bear = solve(dB,B)
    # remove from being in array form
    Be = Bear[1]
    # now sub this into Se and We
    Se = SymPy.subs(Se,B=>Be)
    We = SymPy.subs(We,B=>Be)
    # equation for dA/dt
    dA = (k*S*r)/(r + f*B^2) - kmin*A - K*A + Kmin*W
    # sub in B, S and W expressions
    dA = SymPy.subs(dA,B=>Be,S=>Se,W=>We)
    # then sub in parameters
    dA = SymPy.subs(dA,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    dA = SymPy.subs(dA,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
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
    Be = SymPy.subs(Be,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Be = SymPy.subs(Be,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    Se = SymPy.subs(Se,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    Se = SymPy.subs(Se,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    We = SymPy.subs(We,K=>ps[1],k=>ps[2],Q=>ps[3],q=>ps[4],kmin=>ps[5],Kmin=>ps[6],qmin=>ps[7])
    We = SymPy.subs(We,Qmin=>ps[8],f=>ps[9],r=>ps[10],F=>ps[11],N=>ps[12])
    # then sub the various As values to get final values
    Bs = zeros(3)
    Ss = zeros(3)
    Ws = zeros(3)
    for i = 1:3
        Bs[i] = SymPy.subs(Be,A=>As[i]) |> float
        Ss[i] = SymPy.subs(Se,A=>As[i]) |> float
        Ws[i] = SymPy.subs(We,A=>As[i]) |> float
    end
    # finally put into states
    ss1 = [As[1],Bs[1],Ss[1],Ws[1]]
    sad = [As[2],Bs[2],Ss[2],Ws[2]]
    ss2 = [As[3],Bs[3],Ss[3],Ws[3]]
    println(ss1)
    println(sad)
    println(ss2)
    flush(stdout)
    return(ss1,sad,ss2)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::Array{Float64,2},NG::Int64,Nmid::Int64)
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

# function to generate the variables needed for a given algoritm iteration
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, Ne ]
function genvars(x::Array{Float64,2},λ::SymEngine.Basic,ϑ::Array{SymEngine.Basic,1},NG::Int64,Nmid::Int64)
    # define neccesary symbols
    A, B, S, W, y1, y2, y3, y4 = SymEngine.symbols("A B S W y1 y2 y3 y4")
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
        λt = SymEngine.subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
        λs[i] = SymEngine.subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
    end
    λs[Nmid] = 0 # midpoint should be a saddle point
    for i = Nmid+1:NG
        λt = λ # temporary λ to avoid changing the master one
        λt = SymEngine.subs(λt, A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
        λs[i] = SymEngine.subs(λt, y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
    end
    # Critical points so expect both to be zero
    λs[1] = 0
    λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 4)
    ϑt = Array{SymEngine.Basic,1}(undef,4)
    for j = 1:4
        for i = 2:NG
            ϑt[j] = SymEngine.subs(ϑ[j], A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
            ϑs[i,j] = SymEngine.subs(ϑt[j], y1=>xprim[i,1], y2=>xprim[i,2], y3=>xprim[i,3], y4=>xprim[i,4]) |> float
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
function g!(F::Array{Float64,2},x::Array{Float64,2},C::Array{Float64,1},K::Array{Float64,2},xi::Array{Float64,2},NG::Int64,Nmid::Int64)
    for j = 1:4
        # Start point
        F[1,j] = x[1,j] - xi[1,j]
        # Mid point
        F[Nmid,j] = x[Nmid,j] - xi[2,j]
        # end point
        F[NG+1,j] = x[NG+1,j] - xi[3,j]
    end

    # first path
    for i = 2:Nmid-1
        for j = 1:4
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    # second path
    for i = Nmid+1:NG
        for j = 1:4
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    return(F)
end

# function to solve the system of linear equations
function linsys(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},
                Hx::Array{SymEngine.Basic,1},Hθ::Array{SymEngine.Basic,1},Hθθ::Array{SymEngine.Basic,2},
                Hθx::Array{SymEngine.Basic,2},Δτ::Float64,NG::Int64,Nmid::Int64,H::SymEngine.Basic)
    # define relevant symbols
    A, B, S, W, the1, the2, the3, the4 = SymEngine.symbols("A B S W the1 the2 the3 the4")
    # Make array to store fixed points
    xi = fill(NaN, 3, 4)
    # the fixed points are allowed to vary as both are at zeros
    # Start point
    Hθt = Array{SymEngine.Basic,1}(undef,4) # temporary hamiltonian so master isn't changed
    Hxt = Array{SymEngine.Basic,1}(undef,4)
    for i = 1:4
        Hθt[i] = SymEngine.subs(Hθ[i], the1=>0.0, the2=>0.0, the3=>0.0, the4=>0.0)
        Hxt[i] = SymEngine.subs(Hx[i], the1=>0.0, the2=>0.0, the3=>0.0, the4=>0.0)
    end
    Hθtt = Array{SymEngine.Basic,1}(undef,4)
    Hθttt = Array{SymEngine.Basic,1}(undef,4)
    Ht = SymEngine.subs(H, the1=>0.0, the2=>0.0, the3=>0.0, the4=>0.0)
    Ht = SymEngine.subs(Ht, A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3], W=>x[Nmid,4]) |> float
    # loop to define fixed points
    for i = 1:4
        Hxt[i] = SymEngine.subs(Hxt[i], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3], W=>x[Nmid,4]) |> float
        Hθttt[i] = SymEngine.subs(Hθt[i], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3], W=>x[Nmid,4]) |> float
        Hθtt[i] = SymEngine.subs(Hθt[i], A=>x[1,1], B=>x[1,2], S=>x[1,3], W=>x[1,4]) |> float
        Hθt[i] = SymEngine.subs(Hθt[i], A=>x[end,1], B=>x[end,2], S=>x[end,3], W=>x[end,4]) |> float
        # start point
        xi[1,i] = Δτ*(Hθtt[i]) + x[1,i]
        # midpoint
        xi[2,i] = -Δτ*(Hxt[i]*Ht) + x[Nmid,i]
        # End point
        xi[3,i] = Δτ*(Hθt[i]) + x[end,i]
    end
    Hxθt = Array{SymEngine.Basic,2}(undef,4,4)
    # loop to calculate Hxθ
    for j = 1:4
        for i = 1:4
            Hxθt[i,j] = SymEngine.subs(Hθx[i,j], the1=>0.0, the2=>0.0, the3=>0.0, the4=>0.0)
            Hxθt[i,j] = SymEngine.subs(Hxθt[i,j], A=>x[Nmid,1], B=>x[Nmid,2], S=>x[Nmid,3], W=>x[Nmid,4]) |> float
        end
    end
    # now transpose to get right form
    Hxθt = transpose(Hxθt) # possible source of error if ive misunderstood here
    # loop for additional midpoint terms
    for j = 1:4
        for i = 1:4
            xi[2,i] -= Δτ*(Hxθt[i,j]*Hθttt[j])
        end
    end
    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 4)
    Hxt = Array{SymEngine.Basic,1}(undef,4)
    Hθθt = Array{SymEngine.Basic,2}(undef,4,4)
    Hθxt = Array{SymEngine.Basic,2}(undef,4,4)
    # Put initial values for K in
    for j = 1:4
        for i = 2:NG
            K[i,j] = x[i,j]
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end
    for l = 1:4
        for i = 2:NG
            # Save temporary Hamiltonians so that the substitution doesn't overwrite orginal
            Hxt[l] = SymEngine.subs(Hx[l], A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
            Hxt[l] = SymEngine.subs(Hxt[l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
            for m = 1:4
                Hθθt[m,l] = SymEngine.subs(Hθθ[m,l], A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
                Hθθt[m,l] = SymEngine.subs(Hθθt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
                Hθxt[m,l] = SymEngine.subs(Hθx[m,l], A=>x[i,1], B=>x[i,2], S=>x[i,3], W=>x[i,4])
                Hθxt[m,l] = SymEngine.subs(Hθxt[m,l], the1=>ϑs[i,1], the2=>ϑs[i,2], the3=>ϑs[i,3], the4=>ϑs[i,4]) |> float
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

# main function generates symbolic matrices that the 4 variables can be subbed into
# make vector of these optimisation parameters
function gMAP(ps::Array{Float64,1},NM::Int64,NG::Int64,Nmid::Int64,Δτ::Float64,high2low::Bool)
    # confusingly paras is different to ps, this is something I should probably fix but can't be bothered currently
    # generate symbolic forms for equations required for the simulation
    ϑ, λ, Hθ, Hθx, Hθθ, Hx, H = gensyms(ps)
    # Now generate an initial path to optimize over
    ss1, sad, ss2 = nullcline(ps,15.0,0.1)
    a1 = collect(range(ss1[1],stop=sad[1],length=Nmid))
    a2 = collect(range(sad[1],stop=ss2[1],length=NG+2-Nmid))
    a = vcat(a1,a2[2:length(a2)])
    b1 = collect(range(ss1[2],stop=sad[2],length=Nmid))
    b2 = collect(range(sad[2],stop=ss2[2],length=NG+2-Nmid))
    b = vcat(b1,b2[2:length(b2)])
    s1 = collect(range(ss1[3],stop=sad[3],length=Nmid))
    s2 = collect(range(sad[3],stop=ss2[3],length=NG+2-Nmid))
    s = vcat(s1,s2[2:length(s2)])
    w1 = collect(range(ss1[4],stop=sad[4],length=Nmid))
    w2 = collect(range(sad[4],stop=ss2[4],length=NG+2-Nmid))
    w = vcat(w1,w2[2:length(w2)])
    x = hcat(a,b,s,w)
    # reverse direction of x if high2low is false
    if high2low == false
        x = x[end:-1:1,:]
    end
    # Then appropriatly discretise the path such that it works with this algorithm
    x = discretise(x,NG,Nmid)
    # Set up method to tell if is converged
    convrg = false
    l = 0
    xold = x
    while convrg == false
        x, xprim, λs, ϑs, λprim = genvars(x,λ,ϑ,NG,Nmid)
        newx = linsys(x,xprim,λs,ϑs,λprim,Hx,Hθ,Hθθ,Hθx,Δτ,NG,Nmid,H)
        xn = discretise(newx.zero,NG,Nmid)
        S = Ŝ(xn,xprim,λs,ϑs,λprim,NG)
        δ = 0
        for i = 1:NG+1
            for j = 1:4
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        t2 = time_ns()/10^9
        println("$(sum(S)),$(δ)")
        flush(stdout) # needed to get output in log file
        if l % 100 == 0
            plot(x[:,1],x[:,2])
            savefig("../Results/GraphAB$(l).png")
            plot(x[:,3],x[:,4])
            savefig("../Results/GraphSW$(l).png")
            plot(S)
            savefig("../Results/S$(l).png")
            plot(ϑs)
            savefig("../Results/vartheta$(l).png")
            plot(λs)
            savefig("../Results/lambdas$(l).png")
            plot(x[:,1].+x[:,2].+x[:,3].+x[:,4])
            savefig("../Results/total$(l).png")
        end
        # Now overwrite old x
        x = xn
        if l == 1000
            convrg = true
        end
        l += 1
    end
    return(x)
end

# Function to calculate the action of a given path
function Ŝ(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},NG::Int64)
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
function times(x::Array{Float64,2},xprim::Array{Float64,2},λs::Array{Float64,1},ϑs::Array{Float64,2},λprim::Array{Float64,1},NG::Int64)
    ts = zeros(NG+1)
    for i = 2:NG+1
        ts[i] = ts[i-1] + (1/(2*λs[i-1]) + 1/(2*λs[i]))/NG
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::Array{Float64,1},x::Array{Float64,2},NM::Int64)
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

function main()
    # General parameters
    Ω = 1 # system size, this is just a fudge to get my Euler-Maruyama algorithm (later) to work
    K = 1.0
    k = 1.0
    Q = 1.0
    q = 11.0/15
    kmin = 0.5 # now reverse creation is an important process
    qmin = 0.1
    f = 1.0/(Ω^2) # Promoter switching
    r = 10.0
    F = 10.0*Ω
    Kmin = 10.0^-10 # remains neligable though
    Qmin = 10.0^-10
    Ne = 150.0*Ω # number of elements in the system
    # make vector to write out
    ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, F, Ne ]
    # Optimisation parameters
    NM = 600 # number of segments to discretise MAP onto
    NG = 600 # number of segments to optimize gMAP over
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.1 # I've made this choice arbitarily, too large and the algorithm breaks
    high2low = true # Set if starting from high state or low state
    # Now call simulation function with these parameters
    path = gMAP(ps,NM,NG,Nmid,Δτ,high2low)
    plot(path[:,1],path[:,2])
    savefig("../Results/Graph1.png")
    plot(path[:,3],path[:,4])
    savefig("../Results/Graph2.png")
    path2 = gMAP(ps,NM,NG,Nmid,Δτ,~high2low)
    plot(path[:,1],path[:,2])
    savefig("../Results/Graph12.png")
    plot(path[:,3],path[:,4])
    savefig("../Results/Graph22.png")
    # Now print out the path
    # Block of code to write all this data to a file so I can go through it
    if length(ARGS) >= 1
        output_file = "../Results/$(ARGS[1])1.csv"
        out_file = open(output_file, "w")
        # open file for writing
        for i = 1:size(path,1)
            line = "$(path[i,1]),$(path[i,2]),$(path[i,3]),$(path[i,4])\n"
            write(out_file, line)
        end
        # then close file
        close(out_file)
        output_file2 = "../Results/$(ARGS[1])2.csv"
        out_file = open(output_file2, "w")
        # open file for writing
        for i = 1:size(path2,1)
            line = "$(path2[i,1]),$(path2[i,2]),$(path2[i,3]),$(path2[i,4])\n"
            write(out_file, line)
        end
        # then close file
        close(out_file)
        output_filep = "../Results/$(ARGS[1])p.csv"
        out_filep = open(output_filep, "w")
        for i = 1:length(ps)
            line = "$(ps[i])\n"
            write(out_filep,line)
        end
        # then close file
        close(out_filep)
    end
end


@time main()
