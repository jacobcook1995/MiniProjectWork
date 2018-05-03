#!/usr/bin/env julia
# fr.jl
# A script to plot the unbound ratio for a given f and r over A

# Add package for plotting
using Plots
import GR

# Parameters
const Ω = 50 # system size
const K = 10
const k = K*Ω # steady state for A=k/K=1
const Q = 1
const q = Q*Ω
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const Kmin = 10.0^-20
const qmin = 10.0^-20
const Qmin = 10.0^-20
const f = 25000/(Ω^2) # Promoter switching
const r = 10

const a = 2
const b = 2
const LA = Ω # side length of A axis of the grid
const LB = Ω # side length of B axis of the grid
const δ = 0.1 # detail of grid
const δ2 = 0.01 # detail of nullclines

# function to calculate switching ratio
function ratio(A, B)
    ratioa = r/(r + f*B*B)
    ratiob = r/(r + f*A*A)
    return(ratioa, ratiob)
end

# Produce nullcline functions to be plotted
function nullclines(B)
    A1 = k*r/K./(r+f*B.^a)
    A2 = (r/f*(q./(Q*B)-1)).^(1/b)
    return(A1, A2)
end

function main()
    # A = 0:0.1:Ω
    # B = Ω:-0.1:0
    # ratio = zeros(length(A),2)
    # ratio[:,1] = r./(r + f*A.*A)
    # ratio[:,2] = r./(r + f*B.*B)
    # plot(ratio)
    # savefig("../Results/Ratio.png")
    # Heat map stuff
    sA = convert(Int64, LA/δ) + 1 # Step only works well if L/δ is an integer
    sB = convert(Int64, LB/δ) + 1
    arr = zeros(sA,sB,2)

    # Matrix of A and B values
    for i = 1:sA
        for j = 1:sB
            arr[i,j,1] = (i - 1)*δ
            arr[i,j,2] = (j - 1)*δ
        end
    end

    # More detailed stuff
    δA = zeros(sA, sB)
    δB = zeros(sA, sB)
    for i = 1:sA
        for j = 1:sB
            rA, rB = ratio(arr[i,j,1], arr[i,j,2])
            δA[i,j] = k*rA - K*arr[i,j,1] # first element A, second element B
            δB[i,j] = q*rB - Q*arr[i,j,2]
        end
    end

    # B2 is used for plotting nullclines
    B2 = 0:δ2:LB

    # find nullclines to plot
    s2 = convert(Int64, LB/δ2) + 1
    A12 = zeros(s2,2)
    A12[:,1], A12[:,2] = nullclines(B2)

    # Set plotting back end
    gr()
    Aticks = collect(linspace(0,LA,sA))
    Bticks = collect(linspace(0,LB,sB))
    # Make heatmap of change in A
    heatmap(Bticks, Aticks, δA, xlabel = "B", ylabel = "A")#, xticks = Bticks, yticks = Aticks)
    plot!(A12, B2, xlim = (0,LA), ylim = (0,LB), legend = false)
    savefig("../Results/HeatMapA.png")
    # Make heatmap of change in B
    heatmap(Bticks, Aticks, δB, xlabel = "B", ylabel = "A")
    plot!(A12, B2, xlim = (0,LA), ylim = (0,LB), legend = false)
    savefig("../Results/HeatMapB.png")
    heatmap(Bticks, Aticks, δA-δB, xlabel = "B", ylabel = "A")
    plot!(A12, B2, xlim = (0,LA), ylim = (0,LB), legend = false)
    savefig("../Results/HeatMapA-B.png")
    heatmap(Bticks, Aticks, δB-δA, xlabel = "B", ylabel = "A")
    plot!(A12, B2, xlim = (0,LA), ylim = (0,LB), legend = false)
    savefig("../Results/HeatMapB-A.png")
end

main()
