#!/usr/bin/env julia
# MAPSchlogl.jl
# A script to optimize the path function of the minimum action path in order
# to determine the minimum path action for the 1D Schlogl model

# Putting the relevant imports in
using Optim
using Plots
using Roots
import GR # Need this to stop world age plotting error?

# Firstly should define constants
const k1 = 0.5 # k_{+1}[A]
const K1 = 3.0 # k_{-1}
const k2 = 1.0 # k_{+2}
const K2 = 1.0 # k_{-2}
const B = 4.0
const V = 200
const high2low = false # Set if starting from high state or low state

# Then set parameters of the optimization
const N = 150 # number of segments optimised over

# Multiplicative Guassian noise matrix
function e!(E, x)
    R = [ k1, K1*x, k2*x*x*x, K2*B*x*x]
    E = sqrt(R[1] + R[2] + R[3] + R[4])
    return E
end

# Now construct the three relevant vectors of equations
function f!(F, x)
    R = [ k1, K1*x, k2*x*x*x, K2*B*x*x]
    F = R[1] - R[2] - (R[3] - R[4]) # check
    return F
end

function fX!(FX, x)
    RX = [ 0, K1, k2*3*x*x, K2*B*2*x ]
    FX = RX[1] - RX[2] - (RX[3] - RX[4])
    return FX
end

# Then construct the necessary matrices
function D!(D, x)
    R = [ k1, K1*x, k2*x*x*x, K2*B*x*x]
    D = sum(R)
    return D
end

function f1!(f1, x)
    R1 = [k1, K1*x]
    f1 = R1[1] - R1[2]
    return f1
end

function f2!(f2, x)
    R2 = [k2*x*x*x, K2*B*x*x]
    f2 = R2[1] - R2[2]
    return f2
end

function D2!(D2, x)
    R = [ k1, K1*x, k2*x*x*x, K2*B*x*x]
    D2 = (sum(R))^2
    return D2
end

function DX!(DA, x)
    RX = [ 0, K1, k2*3*x*x, K2*B*2*x ]
    DX = sum(RX)
    return DX
end

# function to split paths into required form
function Pathdiv(path)
    # split path into vectors a and b
    a = path[:]
    # delete first  and last elements of a and b
    deleteat!(a, 1)
    deleteat!(a, N)
    return(a)
end

# function to recombine paths into a plotable form
function Pathmerge(nonfixed)
    # Should now try to reconstruct the paths
    apath = nonfixed[1:N-1]
    push!(apath,fin) # This should all go in a function if it's repeated
    unshift!(apath,star)
    return(apath)
end

# A function to find the crossing points of the nullclines so they can be used
# as start, end and saddle points
function nullcline()
    # write out equation to be solved
    f(x) = k1 - K1*x - k2*(x^3) + K2*B*(x^2)
    three = false
    Xs = [ 0.0, 0.0, 0.0 ]
    while three == false
        X = fzeros(f, 0, 10, no_pts = 1000)
        if length(X) == 3
            three = true
            Xs = X
        end
    end
    sad = Xs[2]
    if high2low == true
        ss1 = Xs[3]
        ss2 = Xs[1]
    else
        ss1 = Xs[1]
        ss2 = Xs[3]
    end
    println(ss1)
    println(sad)
    println(ss2)
    return(ss1,sad,ss2)
end

function AP(thi, tau) # function to calculate action of a given path
    deltat = tau/N # N segements making up the path, each of equal time
    S = 0 # initialise as zero

    # make an initial d and f
    d = 0.0
    h = 0.0

    for i = 1:N
        if i == 1
            posX = (thi[1] + star)/2
        elseif i == N
            posX = (fin + thi[N-1])/2
        else
            posX = (thi[i] + thi[i-1])/2
        end
        h = f!(h, posX)
        d = D!(d, posX)
        # Here Here Here Here Here Here Here Here Here Here Here Here Here Here Here
        if i == 1
            thiv = (thi[1] - star)/deltat
        elseif i == N
            thiv = (fin - thi[(N-1)])/deltat
        else
            thiv = (thi[i] - thi[i-1])/deltat
        end
        S += (0.5*deltat/d)*((thiv-h)^2)
    end
    return(S)
end

# Should now write a function for the gradient to be used for the iteration
function g!(grads, thi, tau)
    deltat = tau/N
    # probably best to pre-calculate the all the f vectors and D matrices
    fu = zeros(1,N)
    fux = zeros(1,N)
    d = zeros(1,1,N)
    d2 = zeros(1,1,N)
    dx = zeros(1,1,N)
    thidot = zeros(1,N)
    for i = 1:N
        if i == 1
            posX = (thi[1] + star)/2
            thidot[i] = (thi[i] - star)/deltat
        elseif i == N
            posX = (fin + thi[N-1])/2
            thidot[i] = (fin - thi[N-1])/deltat
        else
            posX = (thi[i] + thi[i-1])/2
            thidot[i] = (thi[i] - thi[i-1])/deltat
        end
        fu[:,i] = f!(fu[:,i], posX)
        fux[:,i] = fX!(fux[:,i], posX)
        d[:,:,i] = D!(d[:,:,i], posX)
        d2[:,:,i] = D2!(d2[:,:,i], posX)
        dx[:,:,i] = DX!(dx[:,:,i], posX)
    end

    for i = 1:N-1
        # Do each term on a seperate line for clarity
        grads[i] = 0
        grads[i] += (thidot[i] - fu[i])/d[i]
        grads[i] -= (thidot[i+1] - fu[i+1])/d[i+1]
        grads[i] -= (deltat/2)*(thidot[i+1] - fu[i+1])*(fux[i+1]/d[i+1])
        grads[i] -= (deltat/2)*(thidot[i] - fu[i])*(fux[i]/d[i])
        grads[i] -= (deltat/4)*(dx[i+1]/d2[i+1])*((thidot[i+1]-fu[i+1])^2)
        grads[i] -= (deltat/4)*(dx[i]/d2[i])*((thidot[i]-fu[i])^2)
    end
    return(grads)
end

# function to actually perform the optimisation
function optSt(nonfixed,tau)
    lower = zeros(N-1) # No species can drop below zero
    upper = 1000.0*ones(N-1) # No species can have more than the total amount in the system
    initial_x = nonfixed
    od = OnceDifferentiable(f -> AP(f,tau), (grads, f) -> g!(grads,f,tau), initial_x)
    results = optimize(od, initial_x, lower, upper, Fminbox{LBFGS}(), allow_f_increases = true,
                        iterations = 10000, g_tol = 0.0, f_tol = 0.0, x_tol = 0.0)
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
    α = 4.0
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
        if d <= 10.0^-10 || α < 10.0^-10
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
    ents = zeros(N+1)
    KE = zeros(N+1)
    PE = zeros(N+1)
    acts = zeros(N+1)
    entp = zeros(N+1)
    entf = zeros(N+1)
    h = 0.0
    d = 0.0
    f1 = 0.0
    f2 = 0.0
    deltat = tau/N
    # Remove fixed points from path
    # path = pathmin[2:N]
    for i = 1:N+1 # This gives an extra contribution compared to the optimisation!
        posX = pathmin[i]
        # if i == 1
        #     posX = (star + path[i])/2
        # elseif i == N
        #     posX = (path[i-1] + fin)/2
        # else
        #     posX = (path[i-1] + path[i])/2
        # end
        h = f!(h, posX)
        d = D!(d, posX)
        f1 = f1!(f1, posX)
        f2 = f2!(f2, posX)
        # if i == 1
        #     thiv = (path[i] - star)/deltat
        # elseif i == N
        #     thiv = (fin - path[i-1])/deltat
        # else
        #     thiv = (path[i] - path[i-1])/(deltat)
        # end
        if i == 1
            thiv = (pathmin[i+1] - pathmin[i])/deltat
        elseif i == N + 1
            thiv = (pathmin[i] - pathmin[i-1])/deltat
        else
            # different scheme here
            thiv = (pathmin[i+1] - pathmin[i-1])/(2*deltat)
        end
        ents[i] = h*thiv*deltat/d
        KE[i] = thiv*thiv*deltat/(2*d)
        PE[i] = h*h*deltat/(2*d)
        entp[i] = (f1^2 + f2^2)*deltat/d
        entf[i] = 2*f1*f2*deltat/d
    end
    acts = KE + PE - ents
    return(ents, KE, PE, acts, entp, entf)
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
    times = collect(linspace(0.0,t,N+1))
    tmid = 0
    midi = 0
    for i = 1:N+1
        if pathmin[i] >= mid && high2low == false
            tmid = times[i]
            midi = i
            break
        elseif pathmin[i] <= mid && high2low == true
            tmid = times[i]
            midi = i
            break
        end
    end
    pone = plot(pathmin[:,1], times, lab = "MAP", xaxis = "X", yaxis = "T")
    pone = scatter!(pone, [star], [0.0], lab = "Start", seriescolor = :green) # put start point in
    pone = scatter!(pone, [mid], [tmid], lab = "Saddle", seriescolor = :orange)
    pone = scatter!(pone, [fin], [t], lab = "Finish", seriescolor = :red) # put end point in
    ents, kins, pots, acts, prod, flow =  EntProd(pathmin,t)

    # Block of code to write all this data to a file so I can go through it
    if length(ARGS) >= 1
        output_file = "../Results/$(ARGS[1]).csv"
        out_file = open(output_file, "w")
        # open file for writing
        for i = 1:size(pathmin,1)
            line = "$(pathmin[i,1]),$(pathmin[i,2])\n"
            write(out_file, line)
        end
        # final line written as time and action of MAP
        line = "$(t),$(S)\n"
        write(out_file, line)
        close(out_file)
    end

    act = sum(acts)
    points = [ ents, kins, pots, acts ]
    # segcent = zeros(N,1)
    # for i = 1:N
    #     segcent[i] = (pathmin[i,1] + pathmin[i+1,1])/2
    # end
    ptwo = plot(pathmin, points, xaxis = "A", yaxis = "Entropy Production", marker = :auto)
    ptwo = scatter!(ptwo, [star], [0.0], seriescolor = :green)
    ptwo = scatter!(ptwo, [mid], [0.0], seriescolor = :orange)
    ptwo = scatter!(ptwo, [fin], [0.0], seriescolor = :red, leg = false)
    plot(pone, ptwo, layout = (1,2))
    savefig("../Results/Entropy$(high2low).png")
    points2 = [ ents, prod, flow, (prod-flow)  ]
    pthree = plot(pathmin, points2)
    pthree = scatter!(pthree, [star], [0.0], seriescolor = :green)
    pthree = scatter!(pthree, [mid], [0.0], seriescolor = :orange)
    pthree = scatter!(pthree, [fin], [0.0], seriescolor = :red, leg = false)
    savefig("../Results/$(high2low)Entropies.png")
    println("steady state 1 = $(prod[1]*(N+1)/t)")
    println("saddle = $(prod[midi]*(N+1)/t)")
    println("steady state 2 = $(prod[end]*(N+1)/t)")
end

# Now define the paths
start, sad, finish = nullcline()
const star = start # These are then set constant to allow better optimisation
const mid = sad
const fin = finish

# find saddle points
const thi1 = collect(linspace(star,fin,N+1))
if high2low == true
    @time run(11.43411894948466,5)
else
    @time run(17.05863066823466,5)
end
# From high X state: 11.43411894948466
# From low X state: 17.05863066823466
