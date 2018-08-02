#!/usr/bin/env julia
# readparas.jl
# A script to read in files of forward and backward paths, calculate the path entropy productions
# also should read in parameters in order to construct Shannon entropy, Schnakenberg entropy production
# reaction affinities and fluxes
using Plots
using Roots
import GR

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function nullcline(ps::Array{Float64,1},high2low::Bool)
    # gonna use SymPy here
    # A1(x) = sqrt.((r/f)*(q./((qmin+Q).*x-Qmin) - 1))
    # A2(x) = (1/(kmin+K))*((k*r)./(r+f*x.^2) + Kmin)
    # A2 fine A1 definetly wrong somehow
    A1(x) = real(sqrt(complex((ps[10]/ps[9])*(ps[4]/((ps[3]+ps[7])*x - ps[8]) - 1))))
    A2(x) = (1/(ps[5]+ps[1]))*((ps[2]*ps[10])/(ps[10]+ps[9]*x^2) + ps[6])
    g(x) = A1(x) - A2(x)
    three = false
    n = 0
    bs = []
    while three == false
        bs = fzeros(g, 0, 15.0)
        n = length(bs)
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
    return (ss1,sad,ss2)
end

# Vector of functions from MAP case
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function f!(F::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    F[1] = ps[2]*ps[10]/(ps[10] + ps[9]*x[2]^2) - ps[5]*x[1] - ps[1]*x[1] + ps[6]
    F[2] = ps[4]*ps[10]/(ps[10] + ps[9]*x[1]^2) - ps[7]*x[2] - ps[3]*x[2] + ps[8]
    return F
end

# Diffusion matrix from MAP case
function D!(D::Array{Float64,2},x::Array{Float64,1},ps::Array{Float64,1})
    D[1,1] = ps[2]*ps[10]/(ps[10] + ps[9]*x[2]^2) + ps[5]*x[1] + ps[1]*x[1] + ps[6]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = ps[4]*ps[10]/(ps[10] + ps[9]*x[1]^2) + ps[7]*x[2] + ps[3]*x[2] + ps[8]
    return D
end


# Function to calculate the action of the discretised path in the MAP formulation
function EntProd(pathmin::Array{Float64,2},tau::Float64,NM::Int64,ps::Array{Float64,1})
    # probably easiest to calculate the entropy production at each point in the path
    ents = zeros(NM, 2)
    entsf = zeros(NM, 2)
    entsp = zeros(NM, 2)
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
        h = f!(h,[posA,posB],ps)
        d = D!(d,[posA,posB],ps)
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
    return(ents,KE,PE,acts)
end

# function to compute the shannon entropy at a fixed point
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function shannon(point::Array{Float64,1},ps::Array{Float64,1})
    pa = ps[2]*ps[10]/(ps[10]+ps[9]*(point[2])^2)
    pamin = ps[5]*point[1]
    pb = ps[4]*ps[10]/(ps[10]+ps[9]*(point[1])^2)
    pbmin = ps[7]*point[2]
    da = ps[1]*point[1]
    damin = ps[6]
    db = ps[3]*point[2]
    dbmin = ps[8]
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
    S = S1 + S2 + S3 + S4
    return(S,F1,F2,F3,F4,A1,A2,A3,A4)
end

function main()
    for i = 1:length(ARGS)
        # Assign the first command line argument to a variable called input_file
        input_file1 = "../Results/0208/$(ARGS[i])1.csv"
        input_file2 = "../Results/0208/$(ARGS[i])2.csv"
        input_filep = "../Results/0208/$(ARGS[i])p.csv"
        points1 = Array{Float64}(0,2)
        points2 = Array{Float64}(0,2)
        ps = Array{Float64,1}(0)

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
                A = parse(Float64, line[1:(comma - 1)])
                B = parse(Float64, line[(comma + 1):L])
                points1 = vcat(points1, [ A B ])
            end
        end
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
                A = parse(Float64, line[1:(comma - 1)])
                B = parse(Float64, line[(comma + 1):L])
                points2 = vcat(points2, [ A B ])
            end
        end
        # now do second file
        open(input_filep, "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                p = parse(Float64, line)
                ps = vcat(ps, p)
            end
        end
        t1 = ps[11]
        t2 = ps[12]
        N1 = 600
        N2 = 600
        segcent1 = zeros(N1)
        segcent2 = zeros(N2)
        ents1, kins1, pots1, acts1 = EntProd(points1,t1,N1,ps)
        ents2, kins2, pots2, acts2 = EntProd(points2,t2,N2,ps)
        println("Action1 = $(sum(acts1))")
        println("Action2 = $(sum(acts2))")
        # Got actions so can infer stability, now want the steady state entropy productions
        ss1, sad, ss2 = nullcline(ps,false)
        S1, F11, F12, F13, F14, A11, A12, A13, A14 = shannon(ss1,ps)
        println("EntProd1 = $(S1)")
        #println("Fluxes1 = ($(F11),$(F12),$(F13),$(F14))")
        #println("Affinities1 = ($(A11),$(A12),$(A13),$(A14))")
        S2, F21, F22, F23, F24, A21, A22, A23, A24 = shannon(ss2,ps)
        println("EntProd2 = $(S2)")
        #println("Fluxes2 = ($(F21),$(F22),$(F23),$(F24))")
        #println("Affinities2 = ($(A21),$(A22),$(A23),$(A24))")
    end
    return(nothing)
end

@time main()
