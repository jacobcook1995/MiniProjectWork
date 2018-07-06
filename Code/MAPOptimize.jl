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
const high2low = false # Set if starting from high state or low state

# Then set parameters of the optimization
const N = 150 # number of segments optimised over

# Multiplicative Guassian noise matrix
function e!(E, x)
    E[1,1] = 1#sqrt(k*r/(r+f*x[2]*x[2]) + K*x[1])
    E[1,2] = 0
    E[2,1] = 0
    E[2,2] = 1#sqrt(q*r/(r+f*x[1]*x[1]) + Q*x[2])
    return E
end

# Now construct the three relevant vectors of equations
function f!(F, x)
    F[1] = k*r/(r+f*x[2]*x[2]) - K*x[1]
    F[2] = q*r/(r+f*x[1]*x[1]) - Q*x[2]
    return F
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

function f1!(F1, x)
    F1[1] = k*r/(r+f*x[2]*x[2])
    F1[2] = q*r/(r+f*x[1]*x[1])
    return F1
end
function f1A!(F1A, x)
    F1A[1] = 0
    F1A[2] = -2*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    return F1A
end
function f1B!(F1B, x)
    F1B[1] = -2*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    F1B[2] = 0
    return F1B
end

function f2!(F2, x)
    F2[1] = K*x[1]
    F2[2] = Q*x[2]
    return F2
end
function f2A!(F2A, x)
    F2A[1] = K
    F2A[2] = 0
    return F2A
end
function f2B!(F2B, x)
    F2B[1] = 0
    F2B[2] = Q
    return F2B
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
    unshift!(apath,star[1])
    unshift!(bpath,star[2])
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
    f1 = [0.0; 0.0]
    f2 = [0.0; 0.0]

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
        f1 = f1!(f1, [posA; posB])
        f2 = f2!(f2, [posA; posB])
        for j = 1:2
            if i == 1
                thiv = (thi[1+(j-1)*(N-1)] - star[j])/deltat
            elseif i == N
                thiv = (fin[j] - thi[(N-1)*j])/deltat
            else
                thiv = (thi[i+(j-1)*(N-1)] - thi[i-1+(j-1)*(N-1)])/deltat
            end
            #S += (0.5*deltat/d[j,j])*(thiv - h[j])^2
            S += (0.5*deltat/d[j,j])*(thiv^2 + h[j]^2 - 2*f1[j]^2 - 2*f2[j]^2)
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
    f1u = zeros(2,N)
    f1ua = zeros(2,N)
    f1ub = zeros(2,N)
    f2u = zeros(2,N)
    f2ua = zeros(2,N)
    f2ub = zeros(2,N)
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
        f1u[:,i] = f1!(f1u[:,i], [posA; posB])
        f1ua[:,i] = f1A!(f1ua[:,i], [posA; posB])
        f1ub[:,i] = f1B!(f1ub[:,i], [posA; posB])
        f2u[:,i] = f2!(f2u[:,i], [posA; posB])
        f2ua[:,i] = f2A!(f2ua[:,i], [posA; posB])
        f2ub[:,i] = f2B!(f2ub[:,i], [posA; posB])
        d[:,:,i] = D!(d[:,:,i], [posA; posB])
        d2[:,:,i] = D2!(d2[:,:,i], [posA; posB])
        da[:,:,i] = DA!(da[:,:,i], [posA; posB])
        db[:,:,i] = DB!(db[:,:,i], [posA; posB])
    end

    for i = 1:N-1
        for j = 1:2
            # Do each term on a seperate line for clarity
            grads[i+(N-1)*(j-1)] = 0
            grads[i+(N-1)*(j-1)] += (thidot[j,i])/d[j,j,i]
            grads[i+(N-1)*(j-1)] -= (thidot[j,i+1])/d[j,j,i+1]
            # grads[i+(N-1)*(j-1)] += (thidot[j,i] - fu[j,i])/d[j,j,i]
            # grads[i+(N-1)*(j-1)] -= (thidot[j,i+1] - fu[j,i+1])/d[j,j,i+1]
            # terms are different depending on j, hence this elseif statement
            if j == 1
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[1,i+1])*(fua[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[2,i+1])*(fua[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[1,i])*(fua[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[2,i])*(fua[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[1,i+1])*(f1ua[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[2,i+1])*(f1ua[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[1,i])*(f1ua[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[2,i])*(f1ua[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[1,i+1])*(f2ua[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[2,i+1])*(f2ua[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[1,i])*(f2ua[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[2,i])*(f2ua[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i+1]/d2[1,1,i+1])*(thidot[1,i+1]^2 + fu[1,i+1]^2 - 2*f1u[1,i+1]^2 - 2*f2u[1,i+1]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i+1]/d2[2,2,i+1])*(thidot[2,i+1]^2 + fu[2,i+1]^2 - 2*f1u[2,i+1]^2 - 2*f2u[2,i+1]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i]/d2[1,1,i])*(thidot[1,i]^2 + fu[1,i]^2 - 2*f1u[1,i]^2 - 2*f2u[1,i]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i]/d2[2,2,i])*(thidot[2,i]^2 + fu[2,i]^2 - 2*f1u[2,i]^2 - 2*f2u[2,i]^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fua[1,i+1]/d[1,1,i+1])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fua[2,i+1]/d[2,2,i+1])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fua[1,i]/d[1,1,i])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fua[2,i]/d[2,2,i])
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
            elseif j == 2
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[1,i+1])*(fub[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[2,i+1])*(fub[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[1,i])*(fub[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] += (deltat/2)*(fu[2,i])*(fub[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[1,i+1])*(f1ub[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[2,i+1])*(f1ub[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[1,i])*(f1ub[1,i+1]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f1u[2,i])*(f1ub[2,i+1]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[1,i+1])*(f2ub[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[2,i+1])*(f2ub[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[1,i])*(f2ub[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= deltat*(f2u[2,i])*(f2ub[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i+1]/d2[1,1,i+1])*(thidot[1,i+1]^2 + fu[1,i+1]^2 - 2*f1u[1,i+1]^2 - 2*f2u[1,i+1]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i+1]/d2[2,2,i+1])*(thidot[2,i+1]^2 + fu[2,i+1]^2 - 2*f1u[2,i+1]^2 - 2*f2u[2,i+1]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i]/d2[1,1,i])*(thidot[1,i]^2 + fu[1,i]^2 - 2*f1u[1,i]^2 - 2*f2u[1,i]^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i]/d2[2,2,i])*(thidot[2,i]^2 + fu[2,i]^2 - 2*f1u[2,i]^2 - 2*f2u[2,i]^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fub[1,i+1]/d[1,1,i+1])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fub[2,i+1]/d[2,2,i+1])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fub[1,i]/d[1,1,i])
                # grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fub[2,i]/d[2,2,i])
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                # grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
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
    od = OnceDifferentiable(f -> AP(f,tau), (grads, f) -> g!(grads,f,tau), initial_x)
    results = optimize(od, initial_x, lower, upper, Fminbox{LBFGS}(), allow_f_increases = true,
                        iterations = 10, g_tol = 0.0, f_tol = 0.0, x_tol = 0.0) # 10000
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
    nois1 = zeros(N, 2)
    nois2 = zeros(N, 2)
    nois3 = zeros(N, 2)
    h = [0.0; 0.0]
    h1 = [0.0; 0.0]
    h2 = [0.0; 0.0]
    d = [0.0 0.0; 0.0 0.0]
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
        h1 = f1!(h1, [posA posB])
        h2 = f2!(h2, [posA posB])
        for j = 1:2
            if i == 1
                thiv = (path[i,j] - star[j])/deltat
            elseif i == N
                thiv = (fin[j] - path[i-1,j])/deltat
            else
                thiv = (path[i,j] - path[i-1,j])/(deltat)
            end
            #ents[i,j] = h[j]*thiv*deltat/d[j,j]
            ents[i,j] = (h1[j]^2 + h2[j]^2)*deltat/d[j,j]
            KE[i,j] = thiv*thiv*deltat/(2*d[j,j])
            PE[i,j] = h[j]*h[j]*deltat/(2*d[j,j])
            if j == 1
                nois1[i,j] = deltat*0.5*-K/Ω
                nois2[i,j] = deltat*0.5*(thiv + K*posA - k*r/(r + f*posB^2))*K/(2*Ω*(k*r/(r + f*posB^2) + K*posA))
                nois3[i,j] = deltat*(K^2)/(32*(Ω^2)*(k*r/(r + f*posB^2) + K*posA))
            else
                nois1[i,j] = deltat*0.5*-Q/Ω
                nois2[i,j] = deltat*0.5*(thiv + Q*posB - q*r/(r + f*posA^2))*Q/(2*Ω*(q*r/(r + f*posA^2) + Q*posB))
                nois3[i,j] = deltat*(Q^2)/(32*(Ω^2)*(q*r/(r + f*posA^2) + Q*posB))
            end
        end
    end
    acts = KE + PE - ents# - nois1 + nois2 + nois3
    return(ents, KE, PE, acts, nois1, nois2, nois3)
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
    pone = plot(pathmin[:,1],pathmin[:,2], lab = "MAP", xaxis = "A", yaxis = "B")
    pone = scatter!(pone, [star[1]], [star[2]], lab = "Start", seriescolor = :green) # put start point in
    pone = scatter!(pone, [inflex[1]], [inflex[2]], lab = "Saddle Point", seriescolor = :orange)
    pone = scatter!(pone, [fin[1]], [fin[2]], lab = "Finish", seriescolor = :red) # put end point in
    ents, kins, pots, acts, nois1, nois2, nois3 =  EntProd(pathmin,t)
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
    points = [ (ents[:,1] + ents[:,2]), (kins[:,1] + kins[:,2]), (pots[:,1] + pots[:,2]), (acts[:,1] + acts[:,2])]#,
                #(nois1[:,1] + nois1[:,2]), (nois2[:,1] + nois2[:,2]), (nois3[:,1] + nois3[:,2]) ]
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
    println(sum(ents))
    savefig("../Results/Entropy$(high2low).png")

end

# Now define the paths
start, saddle, finish = nullcline()
const star = start # These are then set constant to allow better optimisation
const inflex = saddle
const fin = finish

# find saddle points
const pa1 = collect(linspace(star[1],inflex[1],(N/2)+1))
const pa2 = collect(linspace(inflex[1],fin[1],(N/2)+1))
const pa = vcat(pa1,pa2[2:length(pa2)])
const pb1 = collect(linspace(star[2],inflex[2],(N/2)+1))
const pb2 = collect(linspace(inflex[2],fin[2],(N/2)+1))
const pb = vcat(pb1,pb2[2:length(pb2)])

# First make a reasonable first guess of the path vector
# const pa = collect(star[1]:((fin[1]-star[1])/N):fin[1])
# const pb = collect(star[2]:((fin[2]-star[2])/N):fin[2])
const thi1 = hcat(pa,pb)

@time run(18.93854129467008,5)#18.952490234375
#@time main()
