#!/usr/bin/env julia
# gMAP4bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to this
# time done for the 4 species case

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
using SymPy
import GR # Need this to stop world age plotting error?

# first should make generic functions that take in a full set of parameters, that
# can perform the machine alegbra and then sub in the parameter set to generate
# numeric values for each of the instances considered

# make a symbolic diffusion matrix
function Ds()
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, f, F = symbols("A,B,S,W,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,f,F")
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
    A, B, S, W, K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r, f, F = symbols("A,B,S,W,K,k,Q,q,kmin,Kmin,qmin,Qmin,f,r,f,F")
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
    the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
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
    A, B, S, W = symbols("A,B,S,W")
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
    the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
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
    the1, the2, the3, the4 = symbols("the1,the2,the3,the4")
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
    A, B, S, W = symbols("A,B,S,W")
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
    y1, y2, y3, y4 = symbols("y1,y2,y3,y4")
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
    λ = sqrt(num/den)
    return(λ)
end

#function to find a symbolic expression for ϑ the arc posistion that ensures a
# hamiltonian value of zero for a given point x
function ϑs()
    # create necessary symbols and and symbolic expressions
    y1, y2, y3, y4 = symbols("y1,y2,y3,y4")
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
# main function generates symbolic matrices that the 4 variables can be subbed into
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
    ϑ = ϑs()
    print("$(ϑ)\n")
    # make vector of these optimisation parameters
    paras = [ K; k; Q; q; kmin; Kmin; qmin; Qmin; f; r; F ]

    # Make symbolic matrices/vectors/equations to sub A,B,S,W into when needed

    # Optimisation parameters
    NM = 150 # number of segments to discretise MAP onto
    NG = 150 # number of segments to optimize gMAP over
    Nmid = convert(Int64, ceil((NG+1)/2))
    Δτ = 0.001 # I've made this choice arbitarily, too large and the algorithm breaks
    high2low = false # Set if starting from high state or low state

end


@time main()
