# MAPOptimize.jl
# A script to optimize the path function of the minimum action path in order
# to determine the minimum path action for the driven bistable gene switch with
# cooperativity.
# This work draws heavily from the work of Ruben Perez-Carrasco et al (2016)

# Putting the relevant imports in
using Optim
using Plots
# Firstly should define constants
const k = 100
const K = 10
const q = 10
const Q = 1
const r = 1000000
const f = 100000 # 100000

# Then set parameters of the optimization
const N = 150 # number of segments optimised over
const star = [0; 10] # start point
const fin = [10; 0] # end point

# find saddle points
const inflex = [2; 2] # point of inflextion
const pa1 = collect(star[1]:(2*(inflex[1]-star[1])/N):inflex[1])
const pa2 = collect(inflex[1]:(2*(fin[1]-inflex[1])/N):fin[1])
const pa = vcat(pa1,pa2[2:length(pa2)])
const pb1 = collect(star[2]:(2*(inflex[2]-star[2])/N):inflex[2])
const pb2 = collect(inflex[2]:(2*(fin[2]-inflex[2])/N):fin[2])
const pb = vcat(pb1,pb2[2:length(pb2)])

# First make a reasonable first guess of the path vector
# const pa = collect(star[1]:((fin[1]-star[1])/N):fin[1])
# const pb = collect(star[2]:((fin[2]-star[2])/N):fin[2])
const thi1 = hcat(pa,pb)

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

function AP(thi, tau) # function to calculate action of a given path
    deltat = tau/N
    S = 0 # initialise as zero

    # make an initial d and f
    d = [0.0 0.0; 0.0 0.0]
    h = [0.0; 0.0]
    # need a step to massively disfavour paths that go negative
    negs = false
    for i = 1:N-1
        for j = 1:2
            if thi[i+(N-1)*(j-1)] < 0
                negs = true
            end
        end
    end
    if negs == true
        S = 1000000000000 # large number to disfavour path
    end

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
            S += (0.5*deltat/d[j,j])*((thiv-h[j])^2)
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
    results = optimize(f -> AP(f,tau), (grads, f) -> g!(grads,f,tau), nonfixed, LBFGS(),
                       Optim.Options(g_tol = 0.0, f_tol = 0.0, x_tol = 0.0,
                       iterations = 10000, allow_f_increases = true))
    # Get results out of optimiser
    result = Optim.minimizer(results)
    S = Optim.minimum(results)
    print(Optim.iterations(results))
    print("\n")
    print(Optim.iteration_limit_reached(results))
    print("\n")
    print(Optim.x_converged(results))
    print("\n")
    print(Optim.f_converged(results))
    print("\n")
    print(Optim.g_converged(results))
    print("\n")
    print(Optim.converged(results))
    print("\n")
    print(Optim.f_calls(results))
    print("\n")
    print(Optim.g_calls(results))
    print("\n")
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
    α = 1
    β = 0.1
    while true
        h = 1/100
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
            print(d)
            print("\n")
            print(t)
            print(",")
            print(f)
            print("\n")
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
            if fa <= f - α*β*d # commented out as it seems to be stopping the cycle from finishing
                t = t1
                print(t)
                print(",")
                print(fa)
                print("\n")
                appstep = true
            else
                α *= 0.5
            end
        end
    end
end


# function to run a full optimization
# takes a starting path for each simulation thi, and an initial guess for tau, h is step size
function run(tau,noit)
    linesear(tau,noit)
end

#@time run(20,5)
taus = [ 25; 50; 75; 100; 125; 150; 175; 200; 225; 250; 275; 300; 325; 350; 375; 400; 425; 450 ]
S = zeros(length(taus),1)
for i = 1:length(taus)
   @time _, S[i] = optSt2(taus[i],5) # more loops don't add as much time as you'd think only really adds time to incorrectly calculated MAPs
end
plot(taus,S)
savefig("../Results/Graph.png")
# result, S = optSt2(17.5,1) # 3.9232519533501384, 3.9232519533494825
# print(S)
# print("\n")
# path = Pathmerge(result)
# plot(path[:,1],path[:,2])
# savefig("../Results/Graph10.png")
# result, S = optSt2(17.5,2) # 23.000229120254517
# print(S)
# print("\n")
# path = Pathmerge(result)
# plot(path[:,1],path[:,2])
# savefig("../Results/Graph20.png")
# result, S = optSt2(17.5,3)
# print(S)
# print("\n")
# path = Pathmerge(result)
# plot(path[:,1],path[:,2])
# savefig("../Results/Graph30.png")
