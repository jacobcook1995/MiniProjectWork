#!/usr/bin/env julia
# read.jl
# A script to read in my file and find distinct points
using Plots
import GR

# Firstly should define constants
const Ω = 300
const K = 10
const k = K*Ω # steady state for A=k/K=1
const Q = 1
const q = Q*Ω
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const Kmin = 10.0^-20
const qmin = 10.0^-20
const Qmin = 10.0^-20
const f = 1000/(Ω^2) # Promoter switching
const r = 10
const high2low = true # Set if starting from high state or low state

# Vector of functions from MAP case
function f!(F, x)
    F[1] = k*r/(r+f*x[2]*x[2]) - K*x[1]
    F[2] = q*r/(r+f*x[1]*x[1]) - Q*x[2]
    return F
end

# Diffusion matrix from MAP case
function D!(D, x)
    D[1,1] = k*r/(r+f*x[2]*x[2]) + K*x[1]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = q*r/(r+f*x[1]*x[1]) + Q*x[2]
    return D
end

# Function to calculate the action of the discretised path in the MAP formulation
function EntProd(pathmin,tau,NM)
    # probably easiest to calculate the entropy production at each point in the path
    ents = zeros(NM, 2)
    KE = zeros(NM, 2)
    PE = zeros(NM, 2)
    acts = zeros(NM, 2)
    h = [0.0; 0.0]
    d = [0.0 0.0; 0.0 0.0]
    deltat = tau/NM
    # Remove fixed points from path
    path = pathmin[2:NM,:]
    for i = 1:NM # This gives an extra contribution compared to the optimisation!
        if i == 1
            posA = (pathmin[1,1] + path[i,1])/2
            posB = (pathmin[1,2] + path[i,2])/2
        elseif i == NM
            posA = (path[i-1,1] + pathmin[NM+1,1])/2
            posB = (path[i-1,2] + pathmin[NM+1,2])/2
        else
            posA = (path[i-1,1] + path[i,1])/2
            posB = (path[i-1,2] + path[i,2])/2
        end
        h = f!(h, [posA posB])
        d = D!(d, [posA posB])
        for j = 1:2
            if i == 1
                thiv = (path[i,j] - pathmin[1,j])/deltat
            elseif i == NM
                thiv = (pathmin[NM+1,j] - path[i-1,j])/deltat
            else
                thiv = (path[i,j] - path[i-1,j])/(deltat)
            end
            ents[i,j] = h[j]*thiv*deltat/d[j,j]
            KE[i,j] = thiv*thiv*deltat/(2*d[j,j])
            PE[i,j] = h[j]*h[j]*deltat/(2*d[j,j])
        end
    end
    acts = KE + PE - ents
    return(ents, KE, PE, acts)
end

function main()
    # Assign the first command line argument to a variable called input_file
    input_file1 = "../Results/$(ARGS[1]).csv"
    input_file2 = "../Results/$(ARGS[2]).csv"
    points1 = Array{Float64}(0,2)
    points2 = Array{Float64}(0,2)
    t1 = 0
    S1 = 0
    t2 = 0
    S2 = 0
    lincount = 1 # line counter
    no_lins = countlines(input_file1)
    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for i = 1:L
                if line[i] == ','
                    comma = i
                end
            end
            # check if is last line
            if lincount == no_lins
                t1 = parse(Float64, line[1:(comma - 1)])
                S1 = parse(Float64, line[(comma + 1):L])
            else
                A = parse(Float64, line[1:(comma - 1)])
                B = parse(Float64, line[(comma + 1):L])
                points1 = vcat(points1, [ A B ])
                lincount += 1
            end
        end
    end
    lincount = 1 # line counter
    no_lins = countlines(input_file2)
    # now do second file
    open(input_file2, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for i = 1:L
                if line[i] == ','
                    comma = i
                end
            end
            # check if is last line
            if lincount == no_lins
                t2 = parse(Float64, line[1:(comma - 1)])
                S2 = parse(Float64, line[(comma + 1):L])
            else
                A = parse(Float64, line[1:(comma - 1)])
                B = parse(Float64, line[(comma + 1):L])
                points2 = vcat(points2, [ A B ])
                lincount += 1
            end
        end
    end
    print("$(t1),$(S1)\n")
    print("$(t2),$(S2)\n")
    N1 = 150
    N2 = 150
    segcent1 = zeros(N1)
    segcent2 = zeros(N2)
    act1 = zeros(N1)
    act2 = zeros(N2)
    p = plot(points1[:,1],points1[:,2],label="MAP")
    p = plot!(p,points2[:,1],points2[:,2],xaxis = "A",yaxis = "B",label="gMAP")
    p = scatter!(p, [points1[1,1]], [points1[1,2]], seriescolor = :green,label="Start")
    p = scatter!(p, [points1[end,1]], [points1[end,2]], seriescolor = :red,label="Finish")
    savefig("../Results/CombinedGraph.png")
    ents1, kins1, pots1, acts1 = EntProd(points1,t1,N1)
    ents2, kins2, pots2, acts2 = EntProd(points2,t2,N2)
    for i = 1:N1
        segcent1[i] = (points1[i,1] + points1[i+1,1])/2
    end
    for i = 1:N2
        segcent2[i] = (points2[i,1] + points2[i+1,1])/2
    end
    act1[:] = acts1[:,1] + acts1[:,2]
    act2[:] = acts2[:,1] + acts2[:,2]
    pone = plot(segcent1, act1, label = "MAP")
    pone = plot!(pone, segcent2, act2, label = "gMAP")
    #annotate!(50, 2, text("MAP Action = $(sum(act1))\n gMAP Action = $(sum(act2))", :left, font(5, "Courier")))
    ptwo = plot(reverse(act1), label = "MAP")
    ptwo = plot!(ptwo, reverse(act2), label = "gMAP")
    plot(pone, ptwo, layout = (1,2))
    savefig("../Results/CombinedGraph2.png")
    print("MAP Action = $(sum(act1))\ngMAP Action = $(sum(act2))\n")
    ts1 = collect(linspace(2*t1,0.5*t1,100))
    A1 = zeros(length(ts1))
    for i = 1:length(ts1)
        _, _, _, act = EntProd(points1,ts1[i],N1)
        A1[i] = sum(act)
    end
    plot(ts1,A1)
    savefig("../Results/Graph1.png")
    ts2 = collect(linspace(2*t2,4*t2,100))
    A2 = zeros(length(ts2))
    for i = 1:length(ts2)
        _, _, _, act = EntProd(points2,ts2[i],N2)
        A2[i] = sum(act)
    end
    plot(ts2,A2)
    savefig("../Results/Graph2.png")
end

@time main()
