#!/usr/bin/env julia
# nullcline_flow.jl
# A julia script to plot the nullclines and the flow, so that the effect of
# changing r and f values on the preffered flow direction can be seen

using Plots # use GR instead of plots so that MeshGrid works
import GR # need this to stop plotting in function bug

const Ω = 10
const k = 100
const K = k/Ω #/s
const q = 10 # asymmetric switches - second switch slower but same steady states.
const Q = q/(Ω) #    As a result, flow diagram very different from symmetric case,
     #    showing increased basin of attraction of fast switch

const f = 1
const r = 10
const a = 2
const b = 2
const LA = Ω # side length of A axis of the grid
const LB = Ω # side length of B axis of the grid
const δ = 0.5 # detail of grid
const δ2 = 0.01 # detail of nullclines

# Find velocity of a point in species space
function velos(A, B)
    Adot = (k*r/(r+f*B^a) - K*A)/100
    Bdot = (q*r/(r+f*A^b) - Q*B)/100
    return(Adot, Bdot)
end

function velosnorm(A, B)
    Adot = k*r/(r+f*B^a) - K*A
    Bdot = q*r/(r+f*A^b) - Q*B
    Adotnorm = Adot/sqrt(Adot*Adot + Bdot*Bdot)
    Bdotnorm = Bdot/sqrt(Adot*Adot + Bdot*Bdot)
    return(Adotnorm, Bdotnorm)
end

# Produce nullcline functions to be plotted
function nullclines(B)
    A1 = k*r/K./(r+f*B.^a);
    A2 = (r/f*(q./(Q*B)-1)).^(1/b);
    return(A1, A2)
end

function main()
    sA = convert(Int64, LA/δ) + 1 # Step only works well if L/δ is an integer
    sB = convert(Int64, LB/δ) + 1
    arr = zeros(sA,sB,2)

    for i = 1:sA
        for j = 1:sB
            arr[i,j,1] = (i - 1)*δ
            arr[i,j,2] = LB - (j - 1)*δ
        end
    end

    # B2 is used for plotting nullclines
    B2 = 0:δ2:LB

    # find nullclines to plot
    s2 = convert(Int64, LB/δ2) + 1
    A12 = zeros(s2,2)
    A12[:,1], A12[:,2] = nullclines(B2)

    gr() # set backend to gr
    plot(A12, B2, xlabel = "A", ylabel = "B", xlim = (0,LA), ylim = (0,LB), legend = false)
    # quiver(A, B, dotA, dotB) # comment out as want to see just directions
    pts = vec(P2[(arr[i,j,1], arr[i,j,2]) for i = 1:sA, j = 1:sB])
    quiver!(pts, quiver = velos)
    savefig("../Results/Graph.png")
end

main() # run main program
