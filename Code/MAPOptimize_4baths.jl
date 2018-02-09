#!/usr/bin/env julia
# MAPOptimize_4baths.jl
# A script to optimize the path function of the minimum action path in order
# to determine the minimum path action for the driven bistable gene switch with
# cooperativity.
# This time it is done for the case of 4 baths, A, B, S, W
# This work draws heavily from the work of Ruben Perez-Carrasco et al (2016)

# Putting the relevant imports in
using Optim
using Plots
using Roots
import GR # Need this to stop world age plotting error?

# Parameters
const Ω = 30 # system size
const ϕ = 0.1 # ratio ϕ = q/k
const k = 100 # steady state for A=k/K=1
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const K = k/Ω # K=k'
const Kmin = 10.0^-20
const q = k*ϕ # steady state for B=q/Q=1
const qmin = 10.0^-20
const Q = q/Ω # Q=q'
const Qmin = 10.0^-20
const f = 10 # Promoter switching
const r = 10
const F = 250 # removal rate
const Ne = 12000 # number of elements in the system

# Then set parameters of the optimization
const N = 150 # number of segments optimised over
const high2low = false # Set if starting from high state or low state

# Now construct the three relevant vectors of equations
# x[1] = A, x[2] = B, x[3] = W, x[4] = S
function f!(F1, x)
    F1[1] = k*x[4]*r/(r + f*x[2]*x[2]) - (K + kmin)*x[1] + Kmin*x[3]
    F1[2] = q*x[4]*r/(r + f*x[1]*x[1]) - (Q + qmin)*x[2] + Qmin*x[3]
    F1[3] = K*(x[1] + ϕ*x[2]) - Kmin*(1 + ϕ)*x[3] - F
    F1[4] = -k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] - q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F
    return F1
end

function fA!(FA, x)
    FA[1] = -(K + kmin)
    FA[2] = -2*x[4]*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    FA[3] = K
    FA[4] = kmin + 2*x[4]*q*r*f*x[1]/((r+f*x[1]*x[1])^2)
    return FA
end

function fB!(FB, x)
    FB[1] = -2*x[4]*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    FB[2] = -(Q + qmin)
    FB[3] = Q
    FB[4] = qmin + 2*x[4]*k*r*f*x[2]/((r+f*x[2]*x[2])^2)
    return FB
end

function fW!(FW, x)
    FW[1] = Kmin
    FW[2] = Qmin
    FW[3] = -(1 + ϕ)*Kmin
    FW[4] = 0
    return FW
end

function fS!(FS, x)
    FS[1] = k*r/(r + f*x[2]*x[2])
    FS[2] = q*r/(r + f*x[1]*x[1])
    FS[3] = 0
    FS[4] = -k*r/(r + f*x[2]*x[2]) - q*r/(r + f*x[1]*x[1])
    return FS
end

# Diffusion matrix containing the noise on each term (squared)
function D!(D, x)
    D[1,1] = k*x[4]*r/(r + f*x[2]*x[2]) + (K + kmin)*x[1] + Kmin*x[3]
    D[2,2] = q*x[4]*r/(r + f*x[1]*x[1]) + (Q + qmin)*x[2] + Qmin*x[3]
    D[3,3] = K*(x[1] + ϕ*x[2]) + Kmin*(1 + ϕ)*x[3] + F
    D[4,4] = k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] + q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F
    D[1,2:4] = D[2,1] = D[2,3:4] = D[3,1:2] = D[3,4] = D[4,1:3] = 0
    return D
end

function D2!(D2, x)
    D2[1,1] = (k*x[4]*r/(r + f*x[2]*x[2]) + (K + kmin)*x[1] + Kmin*x[3])^2
    D2[2,2] = (q*x[4]*r/(r + f*x[1]*x[1]) + (Q + qmin)*x[2] + Qmin*x[3])^2
    D2[3,3] = (K*(x[1] + ϕ*x[2]) + Kmin*(1 + ϕ)*x[3] + F)^2
    D2[4,4] = (k*x[4]*r/(r + f*x[2]*x[2]) + kmin*x[1] + q*x[4]*r/(r + f*x[1]*x[1]) + qmin*x[2] + F)^2
    D2[1,2:4] = D2[2,1] = D2[2,3:4] = D2[3,1:2] = D2[3,4] = D2[4,1:3] = 0
    return D2
end

function DA!(DA, x)
    DA[1,1] = (K + kmin)
    DA[2,2] = -2*q*x[4]*r*x[1]/((r + f*x[1]*x[1])^2)
    DA[3,3] = K
    DA[4,4] = kmin - 2*q*x[4]*r*x[1]/((r + f*x[1]*x[1])^2)
    DA[1,2:4] = DA[2,1] = DA[2,3:4] = DA[3,1:2] = DA[3,4] = DA[4,1:3] = 0
    return DA
end

function DB!(DB, x)
    DB[1,1] = -2*k*x[4]*r*x[2]/((r + f*x[2]+x[2])^2)
    DB[2,2] = (Q + qmin)
    DB[3,3] = Q
    DB[4,4] = qmin - 2*k*x[4]*r*x[2]/((r + f*x[2]*x[2])^2)
    DB[1,2:4] = DB[2,1] = DB[2,3:4] = DB[3,1:2] = DB[3,4] = DB[4,1:3] = 0
    return DB
end

function DW!(DW, x)
    DW[1,1] = Kmin
    DW[2,2] = Qmin
    DW[3,3] = Kmin*(1 + ϕ)
    DW[4,4] = 0
    DW[1,2:4] = DW[2,1] = DW[2,3:4] = DW[3,1:2] = DW[3,4] = DW[4,1:3] = 0
    return DW
end

function DS!(DS, x)
    DS[1,1] = k*r/(r + f*x[2]*x[2])
    DS[2,2] = q*r/(r + f*x[1]*x[1])
    DS[3,3] = 0
    DS[4,4] = k*r/(r + f*x[2]*x[2]) + q*r/(r + f*x[1]*x[1])
    DS[1,2:4] = DS[2,1] = DS[2,3:4] = DS[3,1:2] = DS[3,4] = DS[4,1:3] = 0
    return DS
end

# function to split paths into required form
function Pathdiv(path)
    # split path into vectors a and b
    a = path[:,1]
    b = path[:,2]
    w = path[:,3]
    s = path[:,4]
    # delete first  and last elements of a and b
    deleteat!(a, 1)
    deleteat!(a, N)
    deleteat!(b, 1)
    deleteat!(b, N)
    deleteat!(w, 1)
    deleteat!(w, N)
    deleteat!(s, 1)
    deleteat!(s, N)
    nonfixed = vcat(a,b,w,s)
    return(nonfixed)
end

# function to recombine paths into a plotable form
function Pathmerge(nonfixed)
    # Should now try to reconstruct the paths
    apath = nonfixed[1:N-1]
    bpath = nonfixed[N:2*(N-1)]
    wpath = nonfixed[(2*N - 1):3*(N-1)]
    spath = nonfixed[(3*N - 2):4*(N-1)]
    push!(apath,fin[1]) # This should all go in a function if it's repeated
    push!(bpath,fin[2])
    push!(wpath,fin[3])
    push!(spath,fin[4])
    unshift!(apath,star[1])
    unshift!(bpath,star[2])
    unshift!(wpath,star[3])
    unshift!(spath,star[4])
    path = hcat(apath,bpath,wpath,spath)
    return(path)
end
# funtion to find the zeros of the system should find 3 zeros ss1, ss2 and the saddle point
function Zeros()
    g(x) = F./(1 + ((r + f*x.^2)./(ϕ*r + (f/ϕ)*((F/K - x).^2)))) + K.*x - F
    three = false
    n = 0
    As = []
    while three == false
        As = fzeros(g, 0, 2*F/K, order = 1)
        n = length(As)
        if n == 3
            three = true
        end
    end
    Bs = zeros(n)
    Ss = zeros(n)
    Ws = zeros(n)
    for i = 1:n
        Bs[i] = (1/ϕ)*((F/K) - As[i])
        Ss[i] = F/(k*r*(1/(r + f*Bs[i]^2) + ϕ/(r + f*As[i]^2)))
        Ws[i] = Ne - As[i] - Bs[i] - Ss[i]
    end
    sad = [ As[2]; Bs[2]; Ws[2]; Ss[2] ]
    if high2low == true
        ss1 = [ As[1]; Bs[1]; Ws[1]; Ss[1] ]
        ss2 = [ As[3]; Bs[3]; Ws[3]; Ss[3] ]
    else
        ss1 = [ As[3]; Bs[3]; Ws[3]; Ss[3] ]
        ss2 = [ As[1]; Bs[1]; Ws[1]; Ss[1] ]
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
    deltat = tau/N
    S = 0 # initialise as zero

    # make an initial d and f
    d = [ 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 ]
    h = [ 0.0; 0.0; 0.0; 0.0 ]
    # need a step to massively disfavour paths that go negative
    negs = false
    for i = 1:N-1
        for j = 1:4
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
            posW = (thi[(2*N - 1)] + star[3])/2
            posS = (thi[(3*N - 2)] + star[4])/2
        elseif i == N
            posA = (fin[1] + thi[N-1])/2
            posB = (fin[2] + thi[2*(N-1)])/2
            posW = (fin[3] + thi[3*(N-1)])/2
            posS = (fin[4] + thi[4*(N-1)])/2
        else
            posA = (thi[i] + thi[i-1])/2
            posB = (thi[i + (N-1)] + thi[i + (N-2)])/2
            posW = (thi[i + 2*(N-1)] + thi[i - 1 + 2*(N-1)])/2
            posS = (thi[i + 3*(N-1)] + thi[i - 1 + 3*(N-1)])/2
        end
        h = f!(h, [posA; posB; posW; posS])
        d = D!(d, [posA; posB; posW; posS])
        for j = 1:4
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
    fu = zeros(4,N)
    fua = zeros(4,N)
    fub = zeros(4,N)
    fuw = zeros(4,N)
    fus = zeros(4,N)
    d = zeros(4,4,N)
    d2 = zeros(4,4,N)
    da = zeros(4,4,N)
    db = zeros(4,4,N)
    dw = zeros(4,4,N)
    ds = zeros(4,4,N)
    thidot = zeros(4,N)
    for i = 1:N
        if i == 1
            posA = (thi[1] + star[1])/2
            posB = (thi[N] + star[2])/2
            posW = (thi[(2*N - 1)] + star[3])/2
            posS = (thi[(3*N - 2)] + star[4])/2
            thidot[1,i] = (thi[i] - star[1])/deltat
            thidot[2,i] = (thi[N] - star[2])/deltat
            thidot[3,i] = (thi[(2*N - 1)] - star[3])/deltat
            thidot[4,i] = (thi[(3*N - 2)] - star[4])/deltat
        elseif i == N
            posA = (fin[1] + thi[N-1])/2
            posB = (fin[2] + thi[2*(N-1)])/2
            posW = (fin[3] + thi[3*(N-1)])/2
            posS = (fin[4] + thi[4*(N-1)])/2
            thidot[1,i] = (fin[1] - thi[N-1])/deltat
            thidot[2,i] = (fin[2] - thi[2*(N-1)])/deltat
            thidot[3,i] = (fin[3] - thi[3*(N-1)])/deltat
            thidot[4,i] = (fin[4] - thi[4*(N-1)])/deltat
        else
            posA = (thi[i] + thi[i-1])/2
            posB = (thi[i + (N-1)] + thi[i + (N-2)])/2
            posW = (thi[i + 2*(N-1)] + thi[i - 1 + 2*(N-1)])/2
            posS = (thi[i + 3*(N-1)] + thi[i - 1 + 3*(N-1)])/2
            thidot[1,i] = (thi[i] - thi[i-1])/deltat
            thidot[2,i] = (thi[i + (N-1)] - thi[i + (N-2)])/deltat
            thidot[3,i] = (thi[i + 2*(N-1)] - thi[i - 1 + 2*(N-1)])/deltat
            thidot[4,i] = (thi[i + 3*(N-1)] - thi[i - 1 + 3*(N-1)])/deltat
        end
        fu[:,i] = f!(fu[:,i], [posA; posB; posW; posS])
        fua[:,i] = fA!(fua[:,i], [posA; posB; posW; posS])
        fub[:,i] = fB!(fub[:,i], [posA; posB; posW; posS])
        fuw[:,i] = fW!(fuw[:,i], [posA; posB; posW; posS])
        fus[:,i] = fS!(fus[:,i], [posA; posB; posW; posS])
        d[:,:,i] = D!(d[:,:,i], [posA; posB; posW; posS])
        d2[:,:,i] = D2!(d2[:,:,i], [posA; posB; posW; posS])
        da[:,:,i] = DA!(da[:,:,i], [posA; posB; posW; posS])
        db[:,:,i] = DB!(db[:,:,i], [posA; posB; posW; posS])
        dw[:,:,i] = DW!(dw[:,:,i], [posA; posB; posW; posS])
        ds[:,:,i] = DS!(ds[:,:,i], [posA; posB; posW; posS])
    end

    for i = 1:(N-1)
        for j = 1:4
            # Do each term on a seperate line for clarity
            grads[i+(N-1)*(j-1)] = 0
            grads[i+(N-1)*(j-1)] += (thidot[j,i] - fu[j,i])/d[j,j,i]
            grads[i+(N-1)*(j-1)] -= (thidot[j,i+1] - fu[j,i+1])/d[j,j,i+1]
            # terms are different depending on j, hence this elseif statement
            if j == 1
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fua[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fua[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i+1]-fu[3,i+1])*(fua[3,i+1]/d[3,3,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i+1]-fu[4,i+1])*(fua[4,i+1]/d[4,4,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fua[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fua[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i]-fu[3,i])*(fua[3,i]/d[3,3,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i]-fu[4,i])*(fua[4,i]/d[4,4,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[3,3,i+1]/d2[3,3,i+1])*((thidot[3,i+1]-fu[3,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[4,4,i+1]/d2[4,4,i+1])*((thidot[4,i+1]-fu[4,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[3,3,i]/d2[3,3,i])*((thidot[3,i]-fu[3,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(da[4,4,i]/d2[4,4,i])*((thidot[4,i]-fu[4,i])^2)
            elseif j == 2
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fub[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fub[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i+1]-fu[3,i+1])*(fub[3,i+1]/d[3,3,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i+1]-fu[4,i+1])*(fub[4,i+1]/d[4,4,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fub[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fub[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i]-fu[3,i])*(fub[3,i]/d[3,3,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i]-fu[4,i])*(fub[4,i]/d[4,4,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[3,3,i+1]/d2[3,3,i+1])*((thidot[3,i+1]-fu[3,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[4,4,i+1]/d2[4,4,i+1])*((thidot[4,i+1]-fu[4,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[3,3,i]/d2[3,3,i])*((thidot[3,i]-fu[3,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(db[4,4,i]/d2[4,4,i])*((thidot[4,i]-fu[4,i])^2)
            elseif j == 3
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fuw[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fuw[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i+1]-fu[3,i+1])*(fuw[3,i+1]/d[3,3,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i+1]-fu[4,i+1])*(fuw[4,i+1]/d[4,4,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fuw[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fuw[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i]-fu[3,i])*(fuw[3,i]/d[3,3,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i]-fu[4,i])*(fuw[4,i]/d[4,4,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[3,3,i+1]/d2[3,3,i+1])*((thidot[3,i+1]-fu[3,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[4,4,i+1]/d2[4,4,i+1])*((thidot[4,i+1]-fu[4,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[3,3,i]/d2[3,3,i])*((thidot[3,i]-fu[3,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(dw[4,4,i]/d2[4,4,i])*((thidot[4,i]-fu[4,i])^2)
            elseif j == 4
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i+1]-fu[1,i+1])*(fus[1,i+1]/d[1,1,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i+1]-fu[2,i+1])*(fus[2,i+1]/d[2,2,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i+1]-fu[3,i+1])*(fus[3,i+1]/d[3,3,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i+1]-fu[4,i+1])*(fus[4,i+1]/d[4,4,i+1])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[1,i]-fu[1,i])*(fus[1,i]/d[1,1,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[2,i]-fu[2,i])*(fus[2,i]/d[2,2,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[3,i]-fu[3,i])*(fus[3,i]/d[3,3,i])
                grads[i+(N-1)*(j-1)] -= (deltat/2)*(thidot[4,i]-fu[4,i])*(fus[4,i]/d[4,4,i])
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[1,1,i+1]/d2[1,1,i+1])*((thidot[1,i+1]-fu[1,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[2,2,i+1]/d2[2,2,i+1])*((thidot[2,i+1]-fu[2,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[3,3,i+1]/d2[3,3,i+1])*((thidot[3,i+1]-fu[3,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[4,4,i+1]/d2[4,4,i+1])*((thidot[4,i+1]-fu[4,i+1])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[1,1,i]/d2[1,1,i])*((thidot[1,i]-fu[1,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[2,2,i]/d2[2,2,i])*((thidot[2,i]-fu[2,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[3,3,i]/d2[3,3,i])*((thidot[3,i]-fu[3,i])^2)
                grads[i+(N-1)*(j-1)] -= (deltat/4)*(ds[4,4,i]/d2[4,4,i])*((thidot[4,i]-fu[4,i])^2)
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
    print(results)
    print("\n")
    return(result, S)
end


# This function runs the optimisation multiple times in the hope of reducing the ruggedness of the action landscape
function optSt2(tau,noit)
    notfixed = Pathdiv(thi1)
    result, S = optSt(notfixed, tau) # generate standard path
    for i = 2:noit
        gc() # run garbage collection to free up memory
        result, S = optSt(result, tau)
    end
    pathfin = Pathmerge(result)
    return(pathfin, S)
end

# Now define the paths
start, saddle, finish = Zeros()
const star = start
const inflex = saddle
const fin = finish

# First make a reasonable first guess of the path vector
const pa1 = collect(linspace(star[1],inflex[1],(N/2)+1))
const pa2 = collect(linspace(inflex[1],fin[1],(N/2)+1))
const pa = vcat(pa1,pa2[2:length(pa2)])
const pb1 = collect(linspace(star[2],inflex[2],(N/2)+1))
const pb2 = collect(linspace(inflex[2],fin[2],(N/2)+1))
const pb = vcat(pb1,pb2[2:length(pb2)])
const pw1 = collect(linspace(star[3],inflex[3],(N/2)+1))
const pw2 = collect(linspace(inflex[3],fin[3],(N/2)+1))
const pw = vcat(pw1,pw2[2:length(pw2)])
const ps1 = collect(linspace(star[4],inflex[4],(N/2)+1))
const ps2 = collect(linspace(inflex[4],fin[4],(N/2)+1))
const ps = vcat(ps1,ps2[2:length(ps2)])
const thi1 = hcat(pa,pb,pw,ps)

@time path, S = optSt2(20,5)
print("$(S)\n")
plot(path[:,1],path[:,2])
savefig("../Results/4BathGraph.png")
a =  sum(path,2)
plot(a)
savefig("../Results/4BathGraphTotal.png")
