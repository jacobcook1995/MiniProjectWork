#!/usr/bin/env julia
# readparas2.jl
# A script to read in files of forward and backward paths, calculate the path entropy productions
# also should read in parameters in order to construct Schnakenberg entropy production
# and find forward and backwards actions of path
using Plots
using Roots

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
        bs1 = fzeros(g, 0.0, 0.1) # catches zeros about the origin
        bs2 = fzeros(g, 0.1, 2*ps[4]/ps[3]) # the upper bound here is potentially problematic
        bs = vcat(bs1,bs2)
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

# Vector of functions from MAP case
# ps = [ K, k, Q, q, kmin, Kmin, qmin, Qmin, f, r ]
# ps = [ 1, 2, 3, 4,    5,    6,    7,    8, 9, 10]
function f!(F::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    F[1] = ps[2]*ps[10]/(ps[10] + ps[9]*x[2]^2) - ps[5]*x[1] - ps[1]*x[1] + ps[6]
    F[2] = ps[4]*ps[10]/(ps[10] + ps[9]*x[1]^2) - ps[7]*x[2] - ps[3]*x[2] + ps[8]
    return F
end

function f1!(F::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    F[1] = ps[2]*ps[10]/(ps[10] + ps[9]*x[2]^2) - ps[5]*x[1]
    F[2] = ps[4]*ps[10]/(ps[10] + ps[9]*x[1]^2) - ps[7]*x[2]
    return F
end

function f2!(F::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1})
    F[1] = ps[1]*x[1] - ps[6]
    F[2] = ps[3]*x[2] - ps[8]
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

# ps = [k1, K1, k2, K2, B, V ]
function fS!(F::Float64,x::Float64,ps::Array{Float64,1})
    R = [ ps[1], ps[2]*x, ps[3]*x*x*x, ps[4]*ps[5]*x*x]
    F = R[1] - R[2] - (R[3] - R[4])
    return F
end

# Diffusion matrix from MAP case
function DS!(D::Float64,x::Float64,ps::Array{Float64,1})
    R = [ ps[1], ps[2]*x, ps[3]*x*x*x, ps[4]*ps[5]*x*x]
    D = sum(R)
    return D
end

# Function to calculate the action of the discretised path in the MAP formulation
function EntProd(pathmin::Array{Float64,2},tau::Float64,NM::Int64,ps::Array{Float64,1})
    # probably easiest to calculate the entropy production at each point in the path
    ents = zeros(NM, 2)
    KE = zeros(NM, 2)
    PE = zeros(NM, 2)
    acts = zeros(NM, 2)
    prod = zeros(NM, 2)
    flow = zeros(NM, 2)
    h = [0.0; 0.0]
    h1 = [0.0; 0.0]
    h2 = [0.0; 0.0]
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
        h1 = f1!(h1,[posA,posB],ps)
        h2 = f2!(h2,[posA,posB],ps)
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
            prod[i,j] = (h1[j]^2 + h2[j]^2)*deltat/d[j,j]
            flow[i,j] = 2*h1[j]*h2[j]*deltat/d[j,j]
            # ents[i,j] = h[j]*thiv/d[j,j]
            # KE[i,j] = thiv*thiv/(2*d[j,j])
            # PE[i,j] = h[j]*h[j]/(2*d[j,j])
            # prod[i,j] = (h1[j]^2 + h2[j]^2)/d[j,j]
            # flow[i,j] = 2*h1[j]*h2[j]/d[j,j]
        end
    end
    acts = KE + PE - ents
    return(ents,KE,PE,acts,prod,flow)
end

# Function to calculate the action of the discretised Schlögl path in the MAP formulation
# ps = [k1, K1, k2, K2, B, V ]
function EntProdS(pathmin::Array{Float64,1},tau::Float64,NM::Int64,ps::Array{Float64,1})
    # probably easiest to calculate the entropy production at each point in the path
    ents = zeros(NM)
    KE = zeros(NM)
    PE = zeros(NM)
    acts = zeros(NM)
    h = 0.0
    d = 0.0
    deltat = tau/NM
    # Remove fixed points from path
    path = pathmin[2:NM,:]
    for i = 1:NM # This gives an extra contribution compared to the optimisation!
        if i == 1
            posX = (pathmin[1] + path[i])/2
        elseif i == NM
            posX = (path[i-1] + pathmin[NM+1])/2
        else
            posX = (path[i-1] + path[i])/2
        end
        h = fS!(h,posX,ps)
        d = DS!(d,posX,ps)
        if i == 1
            thiv = (path[i] - pathmin[1])/deltat
        elseif i == NM
            thiv = (pathmin[NM+1] - path[i-1])/deltat
        else
            thiv = (path[i] - path[i-1])/(deltat)
        end
        ents[i] = h*thiv*deltat/d
        KE[i] = thiv*thiv*deltat/(2*d)
        PE[i] = h*h*deltat/(2*d)
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
    S12 = S1 + S2
    S34 = S3 + S4
    S = S12 + S34
    return(S)
end

function graphs1()
    # Assign the first command line argument to a variable called input_file
    input_file1 = "../Results/1809/$(ARGS[1])1.csv"
    input_file2 = "../Results/1809/$(ARGS[1])2.csv"
    input_filep = "../Results/1809/$(ARGS[1])p.csv"
    input_file1S = "../Results/1809/$(ARGS[1])S1.csv"
    input_file2S = "../Results/1809/$(ARGS[1])S2.csv"
    input_filepS = "../Results/1809/$(ARGS[1])Sp.csv"
    points1 = Array{Float64,2}(undef,0,3)
    points2 = Array{Float64,2}(undef,0,3)
    ps = Array{Float64,1}(undef,0)
    # Schlögl stuff as well
    points1S = Array{Float64,2}(undef,0,2)
    points2S = Array{Float64,2}(undef,0,2)
    psS = Array{Float64,1}(undef,0)

    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for j = 1:L
                if line[j] == ','
                    comma = j
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points1 = vcat(points1, [ A B 0.0 ])
        end
    end
    # now do second file
    open(input_file2, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for j = 1:L
                if line[j] == ','
                    comma = j
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points2 = vcat(points2, [ A B 0.0 ])
        end
    end
    # now do third file
    open(input_filep, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            p = parse(Float64, line)
            ps = vcat(ps, p)
        end
    end
    # same for Schlögl
    # Open the input file for reading and close automatically at end
    open(input_file1S, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            X = parse(Float64, line)
            points1S = vcat(points1S, [ X 0.0 ])
        end
    end
    # now do second file
    open(input_file2S, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            X = parse(Float64, line)
            points2S = vcat(points2S, [ X 0.0 ])
        end
    end
    # now do third file
    open(input_filepS, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            p = parse(Float64, line)
            psS = vcat(psS, p)
        end
    end
    t1 = ps[11]
    t2 = ps[12]
    N1 = 600
    N2 = 600
    τ1 = t1/(N1+1)
    τ2 = t2/(N2+1)
    t1S = psS[end-1]
    t2S = psS[end]
    N1S = size(points1S,1) - 1
    N2S = size(points2S,1) - 1
    τ1S = t1S/(N1S+1)
    τ2S = t2S/(N2S+1)
    segs1 = zeros(N1,2)
    segs2 = zeros(N2,2)
    segs1S = zeros(N1S,1)
    segs2S = zeros(N2S,1)
    for i = 1:N1
        segs1[i,1] = (points1[i+1,1] + points1[i,1])/2
        segs1[i,2] = (points1[i+1,2] + points1[i,2])/2
    end
    for i = 1:N2
        segs2[i,1] = (points2[i+1,1] + points2[i,1])/2
        segs2[i,2] = (points2[i+1,2] + points2[i,2])/2
    end
    for i = 1:N1S
        segs1S[i] = (points1S[i+1,1] + points1S[i,1])/2
    end
    for i = 1:N2S
        segs2S[i] = (points2S[i+1,1] + points2S[i,1])/2
    end
    # these formula are wrong
    # Add velocities to array
    points1[1,3] = sqrt(((points1[2,1]-points1[1,1])/(τ1))^2 + ((points1[2,2]-points1[1,2])/(τ1))^2)
    for i = 2:N1
        points1[i,3] = sqrt(((points1[i+1,1]-points1[i-1,1])/(2*τ1))^2 + ((points1[i+1,2]-points1[i-1,2])/(2*τ1))^2)
    end
    points1[N1+1,3] = sqrt(((points1[N1+1,1]-points1[N1,1])/(τ1))^2 + ((points1[N1+1,2]-points1[N1,2])/(τ1))^2)
    points2[1,3] = sqrt(((points2[2,1]-points2[1,1])/(τ2))^2 + ((points2[2,2]-points2[1,2])/(τ2))^2)
    for i = 2:N2
        points2[i,3] = sqrt(((points2[i+1,1]-points2[i-1,1])/(2*τ2))^2 + ((points2[i+1,2]-points2[i-1,2])/(2*τ2))^2)
    end
    points2[N2+1,3] = sqrt(((points2[N2+1,1]-points2[N2,1])/(τ2))^2 + ((points2[N2+1,2]-points2[N2,2])/(τ2))^2)
    # And for Schlogl
    points1S[1,2] = abs((points1S[2,1]-points1S[1,1])/(τ1S))
    for i = 2:N1S
        points1S[i,2] = abs((points1S[i+1,1]-points1S[i-1,1])/(2*τ1S))
    end
    points1S[N1S+1,2] = abs((points1S[N1S+1,1]-points1S[N1S,1])/(τ1S))
    points2S[1,2] = abs((points2S[2,1]-points2S[1,1])/(τ2S))
    for i = 2:N2S
        points2S[i,2] = abs((points2S[i+1,1]-points2S[i-1,1])/(2*τ2S))
    end
    points2S[N2S+1,2] = abs((points2S[N2S+1,1]-points2S[N2S,1])/(τ2S))
    # entropy production calculations
    ents1, kins1, pots1, acts1 = EntProd(points1[:,1:2],t1,N1,ps)
    ents2, kins2, pots2, acts2 = EntProd(points2[:,1:2],t2,N2,ps)
    points3 = points1[N1+1:-1:1,:]
    points4 = points2[N2+1:-1:1,:]
    ents3, kins3, pots3, acts3 = EntProd(points3[:,1:2],t1,N1,ps)
    ents4, kins4, pots4, acts4 = EntProd(points4[:,1:2],t2,N2,ps)
    plot([points1[:,1],points4[:,1]],[points1[:,3],points4[:,3]])
    scatter!([points1[1,1]], [0.0], seriescolor = :green)
    scatter!([points1[end,1]], [0.0], seriescolor = :red)
    savefig("../Results/VelocityA.png")
    plot([points2[:,1],points3[:,1]],[points2[:,3],points3[:,3]])
    scatter!([points2[1,1]], [0.0], seriescolor = :green)
    scatter!([points2[end,1]], [0.0], seriescolor = :red)
    savefig("../Results/VelocityB.png")
    plot([points1[:,1],points4[:,1]],[points1[:,2],points4[:,2]])
    scatter!([points1[1,1]], [points1[1,2]], seriescolor = :green)
    scatter!([points1[end,1]], [points1[end,2]], seriescolor = :red)
    savefig("../Results/PathA.png")
    plot([points2[:,1],points3[:,1]],[points2[:,2],points3[:,2]])
    scatter!([points2[1,1]], [points2[1,2]], seriescolor = :green)
    scatter!([points2[end,1]], [points2[end,2]], seriescolor = :red)
    savefig("../Results/PathB.png")
    segs = [ segs1[:,1], segs2[end:-1:1,1]]
    contrse1 = [ sum(ents1,dims=2), sum(ents4,dims=2)]
    plot(segs, contrse1)
    savefig("../Results/Entropies14.png")
    contrsa1 = [ sum(acts1,dims=2), sum(acts4,dims=2)]
    plot(segs, contrsa1)
    savefig("../Results/Actions14.png")
    segs = [ segs1[:,1], segs1[:,1], segs2[end:-1:1,1], segs2[end:-1:1,1]]
    contrskp1 = [ sum(pots1,dims=2), sum(kins1,dims=2), sum(pots4,dims=2), sum(kins4,dims=2)]
    plot(segs, contrskp1)
    savefig("../Results/KinPots14.png")
    segs = [ segs2[:,1], segs1[end:-1:1,1]]
    contrse1 = [ sum(ents2,dims=2), sum(ents3,dims=2)]
    plot(segs, contrse1)
    savefig("../Results/Entropies23.png")
    contrsa1 = [ sum(acts2,dims=2), sum(acts3,dims=2)]
    plot(segs, contrsa1)
    savefig("../Results/Actions23.png")
    segs = [ segs2[:,1], segs2[:,1], segs1[end:-1:1,1], segs1[end:-1:1,1]]
    contrskp1 = [ sum(pots2,dims=2), sum(kins2,dims=2), sum(pots3,dims=2), sum(kins3,dims=2)]
    plot(segs,contrskp1)
    savefig("../Results/KinPots23.png")
    # now for Schlögl model
    ents1S, kins1S, pots1S, acts1S = EntProdS(points1S[:,1],t1S,N1S,psS)
    ents2S, kins2S, pots2S, acts2S = EntProdS(points2S[:,1],t2S,N2S,psS)
    points3S = points1S[N1S+1:-1:1,:]
    points4S = points2S[N2S+1:-1:1,:]
    ents3S, kins3S, pots3S, acts3S = EntProdS(points3S[:,1],t1S,N1S,psS)
    ents4S, kins4S, pots4S, acts4S = EntProdS(points4S[:,1],t2S,N2S,psS)
    plot([points1S[:,1],points4S[:,1]],[points1S[:,2],points4S[:,2]])
    scatter!([points1S[1,1]], [0.0], seriescolor = :green)
    scatter!([points1S[end,1]], [0.0], seriescolor = :red)
    savefig("../Results/SVelocityA.png")
    plot([points2S[:,1],points3S[:,1]],[points2S[:,2],points3S[:,2]])
    scatter!([points2S[1,1]], [0.0], seriescolor = :green)
    scatter!([points2S[end,1]], [0.0], seriescolor = :red)
    savefig("../Results/SVelocityB.png")
    plot([0:(t1S/N1S):t1S,0:(t2S/N2S):t2S],[points1S[:,1],points4S[:,1]])
    scatter!([0], [points4S[1,1]], seriescolor = :green)
    scatter!([t1S,t2S], [points1S[end,1],points4S[end,1]], seriescolor = :red)
    savefig("../Results/SPathA.png")
    plot([0:(t2S/N2S):t2S,0:(t1S/N1S):t1S],[points2S[:,1],points3S[:,1]])
    scatter!([0], [points2S[1,1]], seriescolor = :green)
    scatter!([t2S,t1S], [points2S[end,1],points3S[end,1]], seriescolor = :red)
    savefig("../Results/SPathB.png")
    segs = [ segs1S, segs2S[end:-1:1]]
    contrse1S = [ ents1S, ents4S]
    plot(segs, contrse1S)
    savefig("../Results/SEntropies14.png")
    contrsa1S = [ acts1S, acts4S]
    plot(segs, contrsa1S)
    savefig("../Results/SActions14.png")
    segs = [ segs1S, segs1S, segs2S[end:-1:1], segs2S[end:-1:1]]
    contrskp1S = [ pots1S, kins1S, pots4S, kins4S]
    plot(segs, contrskp1S)
    savefig("../Results/SKinPots14.png")
    segs = [ segs2S, segs1S[end:-1:1]]
    contrse2S = [ ents2S, ents3S]
    plot(segs, contrse2S)
    savefig("../Results/SEntropies23.png")
    contrsa2S = [ acts2S, acts3S]
    plot(segs, contrsa2S)
    savefig("../Results/SActions23.png")
    segs = [ segs2S, segs2S, segs1S[end:-1:1], segs1S[end:-1:1]]
    contrskp2S = [ pots2S, kins2S, pots3S, kins3S]
    plot(segs, contrskp2S)
    savefig("../Results/SKinPots23.png")
    # Got actions so can infer stability, now want the steady state entropy productions
    ss1, sad, ss2 = nullcline(ps,false)
    S1, S11, S21 = shannon(ss1,ps)
    println("EntProd1 = $(S1),$(S11),$(S21)")
    S2, S12, S22 = shannon(ss2,ps)
    println("EntProd2 = $(S2),$(S12),$(S22)")
    return(nothing)
end

function graphs2()
    input_file1 = "../Results/1809/$(ARGS[1])1.csv"
    input_file2 = "../Results/1809/$(ARGS[1])2.csv"
    input_filep = "../Results/1809/$(ARGS[1])p.csv"
    points1 = Array{Float64,2}(undef,0,2)
    points2 = Array{Float64,2}(undef,0,2)
    ps = Array{Float64,1}(undef,0)
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
            comma = 0
            L = length(line)
            for j = 1:L
                if line[j] == ','
                    comma = j
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
            for j = 1:L
                if line[j] == ','
                    comma = j
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points2 = vcat(points2, [ A B ])
        end
    end
    # make array to store magnitude of D
    L = 2000
    X = 9.25
    dx = X/L
    Dmag = zeros(L,L)
    DA = zeros(L,L)
    DB = zeros(L,L)
    fmag = zeros(L,L)
    fA = zeros(L,L)
    fB = zeros(L,L)
    d = zeros(2,2)
    h = zeros(2)
    for j = 1:L
        for i = 1:L
            h = f!(h,[(i-1)*dx,(j-1)*dx],ps)
            fmag[i,j] = sum(h)
            fA[i,j] = h[1]
            fB[i,j] = h[2]
            d = D!(d,[(i-1)*dx,(j-1)*dx],ps)
            Dmag[i,j] = sum(d)
            DA[i,j] = d[1,1]
            DB[i,j] = d[2,2]
        end
    end
    xs = 0:dx:(X-dx)
    # heatmap(xs,xs,fmag,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/Heatmapf.png")
    # heatmap(xs,xs,fA,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/HeatmapfA.png")
    # heatmap(xs,xs,fB,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/HeatmapfB.png")
    # heatmap(xs,xs,Dmag,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/Heatmap.png")
    # heatmap(xs,xs,DA,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/HeatmapA.png")
    # heatmap(xs,xs,DB,xlabel="B",ylabel="A")
    # plot!(points1[:,2],points1[:,1],label="B to A")
    # plot!(points2[:,2],points2[:,1],label="A to B")
    # savefig("../Results/HeatmapB.png")
    # Got actions so can infer stability, now want the steady state entropy production
    ss1, sad, ss2 = nullcline(ps,false)
    SB = shannon(ss1,ps)
    SA = shannon(ss2,ps)
    # extract times and N's
    t1 = ps[end-1]
    t2 = ps[end]
    N1 = size(points1,1) - 1
    N2 = size(points2,1) - 1
    # find entropies along path
    ents1, kins1, pots1, acts1, prod1, flow1 = EntProd(points1[:,1:2],t1,N1,ps)
    ents2, kins2, pots2, acts2, prod2, flow2 = EntProd(points2[:,1:2],t2,N2,ps)
    # find segcents to plot entropies against
    segcent1 = zeros(N1)
    segcent2 = zeros(N2)
    for i = 1:length(segcent1)
        segcent1[i] = (points1[i,1] + points1[i+1,1])/2
    end
    for i = 1:length(segcent2)
        segcent2[i] = (points2[i,1] + points2[i+1,1])/2
    end
    # now have these and should try and plot as two lines
    plot(segcent1,2*N1*(ents1[:,1]+ents1[:,2])/t1,title="B => A",label="Entropy Production",xlabel="A",ylabel="\\Delta S")
    plot!(segcent1,2*N1*(flow1[:,1]+flow1[:,2])/t1,label="''Entropy Flow''")
    plot!(segcent1,2*N1*(prod1[:,1]+prod1[:,2])/t1,label="''Entropy Production''")
    plot!(segcent1,2*N1*(prod1[:,1]+prod1[:,2]-flow1[:,1]-flow1[:,2])/t1,label="Production-Flow")
    scatter!([segcent1[1]], [0.0], seriescolor = :green, label="")
    scatter!([segcent1[end]], [0.0], seriescolor = :red, label="")
    savefig("../Results/EntsB.png")
    hline!([SA],label="high A Schnakenberg")
    hline!([SB],label="high B Schnakenberg")
    savefig("../Results/EntsBSchnakenberg.png")
    plot(segcent2,2*N2*(ents2[:,1]+ents2[:,2])/t2,title="A => B",label="Entropy Production",xlabel="A",ylabel="\\Delta S")
    plot!(segcent2,2*N2*(flow2[:,1]+flow2[:,2])/t2,label="''Entropy Flow''")
    plot!(segcent2,2*N2*(prod2[:,1]+prod2[:,2])/t2,label="''Entropy Production''")
    plot!(segcent2,2*N2*(prod2[:,1]+prod2[:,2]-flow2[:,1]-flow2[:,2])/t2,label="Production-Flow")
    scatter!([segcent2[1]], [0.0], seriescolor = :green, label="")
    scatter!([segcent2[end]], [0.0], seriescolor = :red, label="")
    savefig("../Results/EntsA.png")
    hline!([SA],label="high A Schnakenberg")
    hline!([SB],label="high B Schnakenberg")
    savefig("../Results/EntsASchnakenberg.png")
    return(nothing)
end

@time graphs2()
