#!/usr/bin/env julia
# gMAP2bath.jl
# A script to compute the geometric minimum action path (gMAP) then want it to
# render the path in the usual (and analysable)

# Putting the relevant imports in
using Plots
using Roots
using NLsolve
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
const NM = 1000 # number of segments to discretise MAP onto
const NG = 150 # number of segments to optimize gMAP over
const Δτ = 0.001 # I've made this choice arbitarily, too large and the algorithm breaks
const Δα = 1/NG

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

# Function for the first differential of the hamiltonian in θ
# x[1,1] = A, x[1,2] = B, x[2,1] = θ₁, x[2,2] = θ₂
function Hθ!(Hθ::AbstractArray, x::AbstractArray)
    Hθ[1] = k*r/(r + f*x[1,2]^2) - K*x[1,1] + x[2,1]*(k*r/(r + f*x[1,2]^2) + K*x[1,1])
    Hθ[2] = q*r/(r + f*x[1,1]^2) - Q*x[1,2] + x[2,2]*(q*r/(r + f*x[1,1]^2) + Q*x[1,2])
    return(Hθ)
end

# Function for the second differential of the hamiltonian in θ
# x[1,1] = A, x[1,2] = B, x[2,1] = θ₁, x[2,2] = θ₂
function Hθθ!(Hθθ::AbstractArray, x::AbstractArray)
    Hθθ[1,1] = k*r/(r + f*x[1,2]^2) + K*x[1,1]
    Hθθ[2,2] = q*r/(r + f*x[1,1]^2) + Q*x[1,2]
    Hθθ[1,2] = Hθθ[2,1] = 0
    return(Hθθ)
end

# Function for the first differential of the hamiltonian in x
# x[1,1] = A, x[1,2] = B, x[2,1] = θ₁, x[2,2] = θ₂
function Hx!(Hx::AbstractArray, x::AbstractArray)
    Hx[1] = K*x[2,1]*((x[2,1]/2) - 1) - q*r*f*x[1,1]*x[2,2]*(2 + x[2,2])/((r + f*x[1,1]^2)^2)
    Hx[2] = Q*x[2,2]*((x[2,2]/2) - 1) - k*r*f*x[1,2]*x[2,1]*(2 + x[2,1])/((r + f*x[1,2]^2)^2)
    return(Hx)
end

# Function for the second differential of the hamiltonian in θ then x
# x[1,1] = A, x[1,2] = B, x[2,1] = θ₁, x[2,2] = θ₂
function Hθx!(Hθx::AbstractArray, x::AbstractArray)
    Hθx[1,1] = K*(x[2,1] - 1)
    Hθx[1,2] = -2*f*x[1,2]*k*r*(1+x[2,1])/((r + f*x[1,2]^2)^2)
    Hθx[2,1] = -2*f*x[1,1]*q*r*(1+x[2,2])/((r + f*x[1,1]^2)^2)
    Hθx[2,2] = Q*(x[2,2] - 1)
    return(Hθx)
end

# function to find λ for a particular vector x and y
# x[1] = A, x[2] = B
function λ(x::AbstractVector,y::AbstractVector)
    # tmps to keep the code somewhat readable
    tmp1 = ((K*x[1] - k*r/(r + f*x[2]^2))^2)*(K*x[1] + k*r/(r + f*x[2]^2)) + ((Q*x[2] - q*r/(r + f*x[1]^2))^2)*(Q*x[2] + q*r/(r + f*x[1]^2))
    tmp2 = (y[1]^2)*(K*x[1] + k*r/(r + f*x[2]^2)) + (y[2]^2)*(Q*x[2] + q*r/(r + f*x[1]^2))
    λ = sqrt(tmp1/tmp2)
    return(λ)
end

# function to find ϑ for a particular vector x, y
# x[1] = A, x[2] = B
function ϑ!(ϑ::AbstractArray,x::AbstractVector,y::AbstractVector,λ=nothing)
    # if a lamda isn't provided by the user calculate one, otherwise skip this step to save computations
    if λ == nothing
        λ = λ(x,y)
    end
    ϑ[1] = (λ*y[1] - k*r/(r + f*x[2]^2) + K*x[1])/(k*r/(r + f*x[2]^2) + K*x[1])
    ϑ[2] = (λ*y[2] - q*r/(r + f*x[1]^2) + Q*x[2])/(q*r/(r + f*x[1]^2) + Q*x[2])
    return(ϑ)
end

# function to generate the variables need to solve the system of linear equations
# x[:,1] = A[:], x[:,2] = B
function genvars(x::AbstractArray)
    # calculate velocities
    xprim = fill(NaN, NG+1, 2)
    for i = 2:NG
        for j = 1:2
            xprim[i,j] = (x[i+1,j] - x[i-1,j])/(2/NG)
        end
    end
    # now find λs
    λs = fill(NaN, NG+1)
    for i = 2:NG
        λs[i] = λ(x[i,:],xprim[i,:])
    end
    # Start and end points are assumed through out code to be critical points
    λs[1] = λs[NG+1] = 0
    # now find ϑs
    ϑs = fill(NaN, NG+1, 2)
    for i = 2:NG
        ϑs[i,:] = ϑ!(ϑs[i,:],x[i,:],xprim[i,:],λs[i])
    end
    # Now find λprim
    λprim = fill(NaN, NG+1)
    for i = 2:NG
        λprim[i] = (λs[i+1] - λs[i-1])/(2/NG)
    end

    return(x,xprim,λs,ϑs,λprim)
end

# function to be solved by NLsolve
function g!(F::AbstractArray, x::AbstractArray, C::AbstractVector, K::AbstractArray, xi::AbstractArray)
    F[1,1] = x[1,1] - xi[1,1]
    F[1,2] = x[1,2] - xi[1,2]
    for i = 2:NG
        for j = 1:2
            F[i,j] = x[i,j] - C[i]*(x[i+1,j] - 2*x[i,j] + x[i-1,j]) - K[i,j]
        end
    end
    F[NG+1,1] = x[NG+1,1] - xi[2,1]
    F[NG+1,2] = x[NG+1,2] - xi[2,2]
    return(F)
end

# function to solve the system of linear equations
function linsys(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector)
    # Make empty arrays/vectors to store differential Hamiltonians
    Hx = zeros(2)
    Hθ = zeros(2)
    Hθθ = zeros(2,2)
    Hθx = zeros(2,2)
    point = zeros(2,2)
    # Make array to store fixed points
    xi = fill(NaN, 2, 2)
    # the fixed points are allowed to vary, this operates as a sanity check on the nucline function
    point[2,:] = [ 0.0; 0.0 ] # zeros for both cases, as at ends of arc
    # Start point
    point[1,:] = x[1,:]
    xi[1,:] = Δτ*(Hθ!(Hθ,point)) + x[1,:]
    # End point
    point[1,:] = x[NG+1,:]
    xi[2,:] = Δτ*(Hθ!(Hθ,point)) + x[NG+1,:]
    # Make vector to store constant terms C
    C = fill(NaN, NG+1)
    for i = 2:NG
        C[i] = Δτ*(λs[i]^2)/(1/(NG^2))
    end
    # Make array to store constant vector K's
    K = fill(NaN, NG+1, 2)
    for i = 2:NG
        # Find differential Hamiltonians at specific point
        point[1,:] = x[i,:]
        point[2,:] = ϑs[i,:]
        Hx = Hx!(Hx,point)
        Hθθ = Hθθ!(Hθθ,point)
        Hθx = Hθx!(Hθx,point)
        for j = 1:2
            K[i,j] = x[i,j]
            # think definition of Hθx here is correct, but could be a source of error in future
            K[i,j] -= Δτ*λs[i]*(Hθx[j,1]*xprim[i,1] + Hθx[j,2]*xprim[i,2])
            # Need to be cleverer when I write in 4 dimensions
            K[i,j] += Δτ*(Hθθ[j,1]*Hx[1] + Hθθ[j,2]*Hx[2])
            K[i,j] += Δτ*λs[i]*λprim[i]*xprim[i,j]
        end
    end

    # Make an initial guess of the path, as prior path
    newxi = x
    # make f! a closure of g! for specific xi, C, K
    f!(F,x) = g!(F,x,C,K,xi)
    # Then put all this into the solver
    newx = nlsolve(f!, newxi)
    return(newx)
end

# function to discretise a path in a way that makes it compatable with the algorithm
function discretise(x::AbstractArray)
    # need to interpolate the data from x onto disx, preserving |x'| = const,
    # i.e equal distance between points
    s = zeros(NG+1)
    s[1] = 0
    for i = 2:NG+1
        dA = x[i,1] - x[i-1,1]
        dB = x[i,2] - x[i-1,2]
        s[i] = s[i-1] + sqrt(dA^2 + dB^2) # Could probably drop the sqrts to speed up the code
    end
    # Divide total arc length into equal segments
    ls = zeros(NG+1)
    for i = 1:NG+1
        ls[i] = (i-1)*s[end]/(NG)
    end
    # Find first index greater than a ls[i] for each i
    inds = fill(0,NG+1)
    j = 1
    for i = 1:NG+1
        higher = false
        while higher == false
            if s[j] >= ls[i] || j == NG+1
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they should be fixed
    disx = zeros(NG+1,2)
    disx[1,:] = x[1,:]
    disx[NG+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NG
        one = inds[i] - 1
        two = inds[i]
        s₀ = s[one]
        s₁ = s[two]
        for j = 1:2
            x₀ = x[one,j]
            x₁ = x[two,j]
            disx[i,j] = x₀ + (ls[i] - s₀)*(x₁ - x₀)/(s₁ - s₀)
        end
    end
    return(disx)
end

# Function to generate an initial path then run the alorithm until a MAP is obtained
# This path is then returned for other functions
function gMAP()
    # First find the steady states and saddle point
    ss1, sad, ss2 = nullcline()
    a1 = collect(linspace(ss1[1],sad[1],(NG/2)+1))
    a2 = collect(linspace(sad[1],ss2[1],(NG/2)+1))
    a = vcat(a1,a2[2:length(a2)])
    b1 = collect(linspace(ss1[2],sad[2],(NG/2)+1))
    b2 = collect(linspace(sad[2],ss2[2],(NG/2)+1))
    b = vcat(b1,b2[2:length(b2)])
    x = hcat(a,b)
    # Then appropriatly discretise the path such that it works with this algorithm
    x = discretise(x)
    # Set up method to tell if is converged
    convrg = false
    l = 0
    while convrg == false
        x, xprim, λs, ϑs, λprim = genvars(x)
        newx = linsys(x,xprim,λs,ϑs,λprim)
        xn = discretise(newx.zero)
        # delta is the sum of the differences of all the points in the path
        δ = 0
        for i = 1:NG+1
            for j = 1:2
                δ += abs(x[i,j] - xn[i,j])
            end
        end
        l += 1
        # Now overwrite old x
        x = xn
        if δ <= 0.000001 # Doesn't work as osciallting between two paths
            # Must be a minimum resolution based on Δτ???
            convrg = true
            print("$(l) steps to converge\n")
        end
    end
    return(x)
end

# Function to calculate the action of a given path
function Ŝ(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector)
    S = 0
    S += (3/(2*NG))*(xprim[2,1]*ϑs[2,1] + xprim[2,2]*ϑs[2,2])
    for i = 3:NG-2
        S += (1/NG)*(xprim[i,1]*ϑs[i,1] + xprim[i,2]*ϑs[i,2])
    end
    S += (3/(2*NG))*(xprim[NG-1,1]*ϑs[NG-1,1] + xprim[NG-1,2]*ϑs[NG-1,2])
    return(S)
end

# function to find the times of each point
function times(x::AbstractArray,xprim::AbstractArray,λs::AbstractVector,ϑs::AbstractArray,λprim::AbstractVector,η)
    # make a trimmed vector of lambdas, i.e times smaller than η are replaced by η
    λt = zeros(NG+1)
    #η = 10.0^-8#10.0^-4
    for i = 1:length(λt)
        if λs[i] >= η
            λt[i] = λs[i]
        else
            λt[i] = η
        end
    end
    ts = zeros(NG+1)
    ts[1] = 0
    for i = 2:NG+1
        ts[i] = ts[i-1] + 1/(2*λt[i-1]) + 1/(2*λt[i])
    end
    return(ts)
end

# function to rediscretise a path from arc discretisation to time discretisation
function timdis(ts::AbstractVector,x::AbstractArray)
    # Make discrete vector of time points
    t = zeros(NM+1)
    for i = 1:NM+1
        t[i] = (i-1)*ts[end]/NM
    end
    # Find index of first element greater than t[i] in ts
    inds = fill(0,NM+1)
    j = 1
    for i = 1:NM+1
        higher = false
        while higher == false
            if ts[j] >= t[i] || j == NG+1
                inds[i] = j
                higher = true
            else
                j += 1
            end
        end
    end
    # First do end points as they are fixed
    path = zeros(NM+1,2)
    path[1,:] = x[1,:]
    path[NM+1,:] = x[NG+1,:]
    # This is done to linear order, which is probably good enough
    for i = 2:NM
        one = inds[i] - 1
        two = inds[i]
        t₀ = ts[one]
        t₁ = ts[two]
        for j = 1:2
            x₀ = x[one,j]
            x₁ = x[two,j]
            path[i,j] = x₀ + (t[i] - t₀)*(x₁ - x₀)/(t₁ - t₀)
        end
    end
    return(path)
end

function main()
    path = gMAP()
    x, xprim, λs, ϑs, λprim = genvars(path)
    # use function Ŝ to find the action associated with this path
    S = Ŝ(x,xprim,λs,ϑs,λprim)
    print("Associated Action = $(S)\n")
    # Now find the times tᵢ
    for i = 2:5
        η = 10.0^-i
        ts = times(x,xprim,λs,ϑs,λprim,η)
        print("Time for full path length = $(ts[end])\n")
        plot(ts)
        savefig("../Results/Graph$(high2low)$(i).png")
        path = timdis(ts,x)
        # Block of code to write all this data to a file so I can go through it
        if length(ARGS) >= 1
            output_file = "../Results/$(ARGS[1])$(i)$(high2low).csv"
            out_file = open(output_file, "w")
            # open file for writing
            for i = 1:size(path,1)
                line = "$(path[i,1]),$(path[i,2])\n"
                write(out_file, line)
            end
            close(out_file)
        end
    end
end

@time main()
