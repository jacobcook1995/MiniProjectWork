#!/usr/bin/env julia
# MAPOptimize.jl
# A script to optimize the path function of the minimum action path in order
# to determine the minimum path action for the driven bistable gene switch with
# cooperativity.
# This work draws heavily from the work of Ruben Perez-Carrasco et al (2016)

# Putting the relevant imports in
using Optim
using Plots
using Roots
import GR # Need this to stop world age plotting error?

# Firstly should define constants
const Ω = 1#300
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
const high2low = false # Set if starting from high state or low state

# Then set parameters of the optimization
const N = 150 # number of segments optimised over

# Multiplicative Guassian noise matrix
function e!(E, x)
    E[1,1] = sqrt(k*r/(r+f*x[2]*x[2]) + K*x[1])
    E[1,2] = 0
    E[2,1] = 0
    E[2,2] = sqrt(q*r/(r+f*x[1]*x[1]) + Q*x[2])
    return E
end

# Now construct the three relevant vectors of equations
function f!(F, x)
    F[1] = k*r/(r+f*x[2]*x[2]) - K*x[1]
    F[2] = q*r/(r+f*x[1]*x[1]) - Q*x[2]
    return F
end

function f1!(f1, x)
    f1[1] = k*r/(r+f*x[2]*x[2]) - kmin
    f1[2] = q*r/(r+f*x[1]*x[1]) - qmin
    return f1
end

function f2!(f2, x)
    f2[1] = K*x[1] - Kmin
    f2[2] = Q*x[2] - Qmin
    return f2
end

function fA!(FA, x)
    FA[1] = -K
    FA[2] = -2*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    return FA
end

function fB!(FB, x)
    FB[1] = -2*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    FB[2] = -Q
    return FB
end

# Then construct the necessary matrices
function D!(D, x)
    D[1,1] = k*r/(r+f*x[2]*x[2]) + K*x[1]
    D[1,2] = 0
    D[2,1] = 0
    D[2,2] = q*r/(r+f*x[1]*x[1]) + Q*x[2]
    return D
end

function D2!(D2, x)
    D2[1,1] = (k*r/(r+f*x[2]*x[2]) + K*x[1])^2
    D2[1,2] = 0
    D2[2,1] = 0
    D2[2,2] = (q*r/(r+f*x[1]*x[1]) + Q*x[2])^2
    return D2
end

function DA!(DA, x)
    DA[1,1] = K
    DA[1,2] = 0
    DA[2,1] = 0
    DA[2,2] = -2*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    return DA
end

function DB!(DB, x)
    DB[1,1] = -2*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    DB[1,2] = 0
    DB[2,1] = 0
    DB[2,2] = Q
    return DB
end

# function to split paths into required form
function Pathdiv(path)
    # split path into vectors a and b
    a = path[:,1]
    b = path[:,2]
    # delete first  and last elements of a and b
    deleteat!(a, 1)
    deleteat!(a, N)
    deleteat!(b, 1)
    deleteat!(b, N)
    nonfixed = vcat(a,b)
    return(nonfixed)
end

# function to recombine paths into a plotable form
function Pathmerge(nonfixed)
    # Should now try to reconstruct the paths
    apath = nonfixed[1:N-1]
    bpath = nonfixed[N:2*(N-1)]
    push!(apath,fin[1]) # This should all go in a function if it's repeated
    push!(bpath,fin[2])
    pushfirst!(apath,star[1])
    pushfirst!(bpath,star[2])
    path = hcat(apath,bpath)
    return(path)
end

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
function nullcline()
    a = 2
    b = 2
    A1(x) = k*r/(K*(r+f*x^a))
    A2(x) = (r/f*(q/(Q*x)-1))^(1/b)
    g(x) = k*r/(K*(r+f*x^a)) - (r/f*(q/(Q*x)-1))^(1/b) #A1(x) - A2(x)
    xs = fzeros(g, 0, q/Q)
    sad = [A1(xs[2]); xs[2]]
    if high2low == true
        ss1 = [A1(xs[1]); xs[1]]
        ss2 = [A1(xs[3]); xs[3]]
    else
        ss1 = [A1(xs[3]); xs[3]]
        ss2 = [A1(xs[1]); xs[1]]
    end
    print(ss1)
    print("\n")
    print(sad)
    print("\n")
    print(ss2)
    print("\n")
    return (ss1,sad,ss2)
end

function AP(thi, tau) # function to calculate action of a given path
    deltat = tau/N # N segements making up the path, each of equal time
    S = 0 # initialise as zero

    # make an initial d and f
    d = [0.0 0.0; 0.0 0.0]
    h = [0.0; 0.0]

    for i = 1:N
        if i == 1
            posA = (thi[1] + star[1])/2
            posB = (thi[N] + star[2])/2
        elseif i == N
            posA = (fin[1] + thi[N-1])/2
            posB = (fin[2] + thi[2*(N-1)])/2
        else
            posA = (thi[i] + thi[i-1])/2
            posB = (thi[i+(N-1)] + thi[i+(N-2)])/2
        end
        h = f!(h, [posA; posB])
        d = D!(d, [posA; posB])
        for j = 1:2
            if i == 1
                thiv = (thi[1+(j-1)*(N-1)] - star[j])/deltat
            elseif i == N
                thiv = (fin[j] - thi[(N-1)*j])/deltat
            else
                thiv = (thi[i+(j-1)*(N-1)] - thi[i-1+(j-1)*(N-1)])/deltat
            end
            S += (0.5*deltat/d[j,j])*(thiv - h[j])^2
        end
    end
    return(S)
end

# Should now write a function for the gradient to be used for the iteration
function g!(grads, thi, tau)
    deltat = tau/N
    # probably best to pre-calculate the all the f vectors and D matrices
    fu = zeros(2,N)
    fua = zeros(2,N)
    fub = zeros(2,N)
    d = zeros(2,2,N)
    d2 = zeros(2,2,N)
    db = zeros(2,2,N)
    da = zeros(2,2,N)
    thidot = zeros(2,N)
    for i = 1:N
        if i == 1
            posA = (thi[1] + star[1])/2
            posB = (thi[N] + star[2])/2
            thidot[1,i] = (thi[i] - star[1])/deltat
            thidot[2,i] = (thi[N] - star[2])/deltat
        elseif i == N
            posA = (fin[1] + thi[N-1])/2
            posB = (fin[2] + thi[2*(N-1)])/2
            thidot[1,i] = (fin[1] - thi[N-1])/deltat
            thidot[2,i] = (fin[2] - thi[2*(N-1)])/deltat
        else
            posA = (thi[i] + thi[i-1])/2
            posB = (thi[i+(N-1)] + thi[i+(N-2)])/2
            thidot[1,i] = (thi[i] - thi[i-1])/deltat
            thidot[2,i] = (thi[i+(N-1)] - thi[i+(N-2)])/deltat
        end
        fu[:,i] = f!(fu[:,i], [posA; posB])
        fua[:,i] = fA!(fua[:,i], [posA; posB])
        fub[:,i] = fB!(fub[:,i], [posA; posB])
        d[:,:,i] = D!(d[:,:,i], [posA; posB])
        d2[:,:,i] = D2!(d2[:,:,i], [posA; posB])
        da[:,:,i] = DA!(da[:,:,i], [posA; posB])
        db[:,:,i] = DB!(db[:,:,i], [posA; posB])
    end

    for i = 1:N-1
        for j = 1:2
            # Do each term on a seperate line for clarity
            grads[i+(N-1)*(j-1)] = 0
            grads[i+(N-1)*(j-1)] += (thidot[j,i] - fu[j,i])/d[j,j,i]
            grads[i+(N-1)*(j-1)] -= (thidot[j,i+1] - fu[j,i+1])/d[j,j,i+1]
            # terms are different depending on j, hence this elseif statement
            if j == 1
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fua[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fua[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fua[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fua[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
            elseif j == 2
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fub[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fub[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fub[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fub[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
            end
        end
    end
    return(grads)
end

# function to actually perform the optimisation
function optSt(nonfixed,tau)
    lower = zeros(2*(N-1)) # No species can drop below zero
    upper = 1000.0*ones(2*(N-1)) # No species can have more than the total amount in the system
    initial_x = nonfixed
    results = optimize(f -> AP(f,tau), (grads, f) -> g!(grads,f,tau), initial_x, LBFGS())
    # od = OnceDifferentiable(f -> AP(f,tau), (grads, f) -> g!(grads,f,tau), initial_x)
    # results = optimize(od, initial_x, lower, upper, Fminbox{LBFGS}(), allow_f_increases = true,
    #                     iterations = 10000, g_tol = 0.0, f_tol = 0.0, x_tol = 0.0)
    # Get results out of optimiser
    result = Optim.minimizer(results)
    S = Optim.minimum(results)
    return(result, S)
end

# This function runs the optimisation multiple times in the hope of reducing the ruggedness of the action landscape
function optSt2(tau,noit)
    notfixed = Pathdiv(thi1)
    result, S = optSt(notfixed, tau) # generate standard path
    for i = 2:noit
        result, S = optSt(result, tau)
    end
    pathfin = Pathmerge(result)
    return(pathfin, S)
end

# A function to run a linesearch
function linesear(tau,noit)
    k = 0 # initialise counter
    t = tau
    α = 4
    β = 0.1
    while true
        h = 1/10000
        # first proper step of the algorithm is to calculate the direction vector
        _, f = optSt2(t,noit)
        _, fh = optSt2(t+h,noit)
        _, fminh = optSt2(t-h,noit)
        dplus = - (fh-f)/h
        dminus = - (fminh-f)/h
        d = max(dplus,dminus)
        if dplus >= dminus
            plus = true
        else
            plus = false
        end
        # Check that the rate of descent is still significant if not end linesearch
        if d <= 10.0^-10
            print("Line Search Done")
            print("\n")
            return(t)
            break
        end
        # Next step is to find the correct step size to make, α,
        appstep = false
        while appstep == false
            if plus == true
                t1 = t + α
            else
                t1 = t - α
            end
            _, fa = optSt2(t1,noit)
            if fa <= f - α*β*d
                t = t1
                appstep = true
                print("Still making progress. τ = $(t)\n")
            else
                α *= 0.5
            end
        end
    end
end

# function to find the entropy production of a path that takes a certain time tau
function EntProd(pathmin,tau)
    # probably easiest to calculate the entropy production at each point in the path
    ents = zeros(N, 2)
    KE = zeros(N, 2)
    PE = zeros(N, 2)
    acts = zeros(N, 2)
    prods = zeros(N, 2)
    flows = zeros(N, 2)
    h = [0.0; 0.0]
    d = [0.0 0.0; 0.0 0.0]
    f1 = [0.0; 0.0]
    f2 = [0.0; 0.0]
    deltat = tau/N
    # Remove fixed points from path
    path = pathmin[2:N,:]
    for i = 1:N # This gives an extra contribution compared to the optimisation!
        if i == 1
            posA = (star[1] + path[i,1])/2
            posB = (star[2] + path[i,2])/2
        elseif i == N
            posA = (path[i-1,1] + fin[1])/2
            posB = (path[i-1,2] + fin[2])/2
        else
            posA = (path[i-1,1] + path[i,1])/2
            posB = (path[i-1,2] + path[i,2])/2
        end
        h = f!(h, [posA posB])
        d = D!(d, [posA posB])
        f1 = f1!(f1, [posA posB])
        f2 = f2!(f2, [posA posB])
        for j = 1:2
            if i == 1
                thiv = (path[i,j] - star[j])/deltat
            elseif i == N
                thiv = (fin[j] - path[i-1,j])/deltat
            else
                thiv = (path[i,j] - path[i-1,j])/(deltat)
            end
            ents[i,j] = h[j]*thiv*deltat/d[j,j]
            KE[i,j] = thiv*thiv*deltat/(2*d[j,j])
            PE[i,j] = h[j]*h[j]*deltat/(2*d[j,j])
            prods[i,j] = (f1[j]^2 + f2[j]^2)*deltat/d[j,j]
            flows[i,j] = 2*f1[j]*f2[j]*deltat/d[j,j]
        end
    end
    acts = KE + PE - ents
    return(ents, KE, PE, acts, prods, flows)
end


# function to run a full optimization
# takes a starting path for each simulation thi, and an initial guess for tau, h is step size
function run(tau,noit)
    t = linesear(tau,noit)
    print(t)
    print("\n")
    pathmin, S = optSt2(t,noit)
    print(S)
    print("\n")
    gr()
    times = collect(range(0.0,stop=t,length=N+1))
    tmid = 0
    midi = 0
    for i = 1:N+1
        if pathmin[i,1] >= inflex[1] && high2low == false
            tmid = times[i]
            midi = i
            break
        elseif pathmin[i,1] <= inflex[1] && high2low == true
            tmid = times[i]
            midi = i
            break
        end
    end
    pone = plot(pathmin[:,1],pathmin[:,2], lab = "MAP", xaxis = "A", yaxis = "B")
    pone = scatter!(pone, [star[1]], [star[2]], lab = "Start", seriescolor = :green) # put start point in
    pone = scatter!(pone, [inflex[1]], [inflex[2]], lab = "Saddle Point", seriescolor = :orange)
    pone = scatter!(pone, [fin[1]], [fin[2]], lab = "Finish", seriescolor = :red) # put end point in
    ents, kins, pots, acts, prod, flow =  EntProd(pathmin,t)
    pathminr = pathmin[end:-1:1,:]
    entsr, kinsr, potsr, actsr, _, _ = EntProd(pathminr,t)
    println("Reverse")
    println(sum(entsr))
    println("Forward")
    println(sum(ents))
    println("Change")
    println(sum(entsr) + sum(ents))
    println(sum(actsr))
    # Block of code to write all this data to a file so I can go through it
    if length(ARGS) >= 1
        output_file = "../Results/$(ARGS[1]).csv"
        out_file = open(output_file, "w")
        # open file for writing
        for i = 1:size(pathmin,1)
            line = "$(pathmin[i,1]),$(pathmin[i,2])\n"
            write(out_file, line)
        end
        # final line written as time and action of gMAP
        line = "$(t),$(S)\n"
        write(out_file, line)
        close(out_file)
    end

    act = sum(acts)
    points = [ (ents[:,1] + ents[:,2]), (kins[:,1] + kins[:,2]), (pots[:,1] + pots[:,2]), (acts[:,1] + acts[:,2])]
    segcent = zeros(N,1)
    for i = 1:N
        segcent[i] = (pathmin[i,1] + pathmin[i+1,1])/2
    end
    ptwo = plot(segcent, points, xaxis = "A", yaxis = "Entropy Production", marker = :auto, legend = false)
    ptwo = scatter!(ptwo, [star[1]], [0.0], seriescolor = :green)
    ptwo = scatter!(ptwo, [inflex[1]], [0.0], seriescolor = :orange)
    ptwo = scatter!(ptwo, [fin[1]], [0.0], seriescolor = :red)#, leg = false)
    #annotate!(inflex[1]+3, minimum(ents)+0.01, text("Action = $(act)\n", :left, font(5, "Courier")))
    plot(pone, ptwo, layout = (1,2))
    println("EntProd = $(sum(ents))")
    println("Change of entropy first half = $(sum(ents[1:midi,1:2]))")
    println("Change of entropy second half = $(sum(ents[(midi+1):end,1:2]))")
    savefig("../Results/Entropy$(high2low).png")
    points2 = [ (ents[:,1] + ents[:,2]), (prod[:,1] + prod[:,2]), (flow[:,1] + flow[:,2]), (prod[:,1] + prod[:,2] - flow[:,1] - flow[:,2])  ]
    pthree = plot(segcent, points2)
    pthree = scatter!(pthree, [star[1]], [0.0], seriescolor = :green)
    pthree = scatter!(pthree, [inflex[1]], [0.0], seriescolor = :orange)
    pthree = scatter!(pthree, [fin[1]], [0.0], seriescolor = :red, leg = false)
    savefig("../Results/$(high2low)Entropies.png")

end

# Now define the paths
start, saddle, finish = nullcline()
const star = start # These are then set constant to allow better optimisation
const inflex = saddle
const fin = finish

# find saddle points
const len = round(Int64,(N/2)+1)
const pa1 = collect(range(star[1],stop=inflex[1],length=len))
const pa2 = collect(range(inflex[1],stop=fin[1],length=len))
const pa = vcat(pa1,pa2[2:length(pa2)])
const pb1 = collect(range(star[2],stop=inflex[2],length=len))
const pb2 = collect(range(inflex[2],stop=fin[2],length=len))
const pb = vcat(pb1,pb2[2:length(pb2)])


const thi1 = hcat(pa,pb)
if high2low == false
    @time run(18.93817579898832,5) # these are valid for N = 150
else
    @time run(21.937873737550667,5)
end
