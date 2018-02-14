#!/usr/bin/env julia
# gillespie_2bath_traj_expl.jl
# Julia Script for Gillespie simulation of toy driven ecosystem, this now
# only calculates the histogram rather than trajectories to save computation time
# In this case 2 resivours are considered, protiens A and B
# This script explicitly generates trajectories so that entropy production can be found

# Add package for plotting
using Plots
import GR

# Parameters
const Ω = 50 # system size, maybe this isn't correct terminology
#const ss = 50 # Steady state size
const K = 10
const k = K*Ω # steady state for A=k/K=1
const Q = 1
const q = Q*Ω
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const Kmin = 10.0^-20
const qmin = 10.0^-20
const Qmin = 10.0^-20
const f = 10/(Ω^2) # Promoter switching
const r = 10

function gillespie()
    Ti = 0.0  # initial time
    Tf = 1000#0 # end of simulation time in s
    batchsize = 100000
    Traj = zeros(0,5)
    Tims = []
    traj = zeros(batchsize,5)
    tims = zeros(batchsize)

    t = Ti # initialise
    A = 0
    B = Ω
    a = 1
    b = 1
    trans = 0

    i = 1
    traj[1,:] = [ A B a b trans ]
    tims[1] = t
    # Main loop
    while t <= Tf
        # set prior A and Bs

        # rates
        rates = [ k*a, kmin*A, K*A, Kmin, r*b, qmin*B, Q*B, Qmin, f*B*(B-1)*a, (1-a)*r, f*A*(A-1)*b, (1-b)*r ]
        rs = rates/sum(rates)
        rone = rand() # first random number

        # which reaction?
        if rone < rs[1]
            A += 1
            trans = 1
        elseif rone < (rs[1]+rs[2])
            A -= 1
            trans = 2
        elseif rone < (rs[1]+rs[2]+rs[3])
            A -= 1
            trans = 3
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4])
            A += 1
            trans = 4
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5])
            B += 1
            trans = 5
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6])
            B -= 1
            trans = 6
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7])
            B -= 1
            trans = 7
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8])
            B += 1
            trans = 8
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8]+rs[9])
            B -= 2
            a -= 1
            trans = 9
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8]+rs[9]+rs[10])
            B += 2
            a += 1
            trans = 10
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8]+rs[9]+rs[10]+rs[11])
            A -= 2
            b -= 1
            trans = 11
        else
            A += 2
            b += 1
            trans = 12
        end

        # update time
        rtwo = rand()  # second random number
        t = t - log(rtwo)/sum(rates)

        i += 1
        traj[i,:] = [ A B a b trans ]
        tims[i] = t
        if i == batchsize
            Traj = vcat(Traj,traj)
            Tims = vcat(Tims,tims)
            traj = zeros(batchsize,5)
            tims = zeros(batchsize)
            i = 0
        end
    end
    # Need to now cat the final incomplete batch
    trajverytemp = fill(0, i, 5)
    timsverytemp = fill(0.0, i)
    # Select the valued entries
    for p = 1:i
        trajverytemp[p,:] = traj[p,:]
        timsverytemp[p] = tims[p]
    end
    clear!(:traj)
    clear!(:tims)
    gc()
    Traj = vcat(Traj, trajverytemp)
    Tims = vcat(Tims, timsverytemp)
    # clear temps
    clear!(:trajverytemp)
    clear!(:timsverytemp)
    gc()

    # Now need to write some code that calculates the entropy production of this trajectory
    L = length(Tims)
    Wf = zeros(L - 1)
    Wminf = zeros(L - 1)
    For = 0
    Back = 0
    for i = 2:L
        # calculate rate of going from pos[i-1] to pos[i]
        rate_type = Traj[i,5]
        if rate_type == 0
            error("You've done something wrong here Jacob.\n")
        elseif rate_type == 1
            Wf[i-1] = k*Traj[i-1,3]
        elseif rate_type == 2
            Wf[i-1] = kmin*Traj[i-1,1]
        elseif rate_type == 3
            Wf[i-1] = K*Traj[i-1,1]
        elseif rate_type == 4
            Wf[i-1] = Kmin
        elseif rate_type == 5
            Wf[i-1] = q*Traj[i-1,4]
        elseif rate_type == 6
            Wf[i-1] = qmin*Traj[i-1,2]
        elseif rate_type == 7
            Wf[i-1] = Q*Traj[i-1,2]
        elseif rate_type == 8
            Wf[i-1] = Qmin
        elseif rate_type == 9
            Wf[i-1] = f*Traj[i-1,2]*(Traj[i-1,2]-1)*Traj[i-1,3]
        elseif rate_type == 10
            Wf[i-1] = (1-Traj[i-1,3])*r
        elseif rate_type == 11
            Wf[i-1] = f*Traj[i-1,1]*(Traj[i-1,1]-1)*Traj[i-1,4]
        elseif rate_type == 12
            Wf[i-1] = (1-Traj[i-1,4])*r
        else
            error("You've done something wrong here Jacob.\n")
        end
        For += log(Wf[i-1])
    end
    # Now find reversed trajectory
    for i = L:(-1):2
        # calculate rate of going from pos[i] to pos[i-1]
        rate_type = Traj[i,5]
        if rate_type == 0
            error("You've done something wrong here Jacob.\n")
        elseif rate_type == 1
            Wminf[i-1] = kmin*Traj[i,1]
        elseif rate_type == 2
            Wminf[i-1] = r*k/(r + f*Traj[i,2]*(Traj[i,2]-1))
        elseif rate_type == 3
            Wminf[i-1] = Kmin
        elseif rate_type == 4
            Wminf[i-1] = K*Traj[i,1]
        elseif rate_type == 5
            Wminf[i-1] = qmin*Traj[i,2]
        elseif rate_type == 6
            Wminf[i-1] = r*q/(r + f*Traj[i,1]*(Traj[i,1]-1))
        elseif rate_type == 7
            Wminf[i-1] = Qmin
        elseif rate_type == 8
            Wminf[i-1] = Q*Traj[i,2]
        elseif rate_type == 9
            Wminf[i-1] = (1-Traj[i,3])*r
        elseif rate_type == 10
            Wminf[i-1] = f*Traj[i,2]*(Traj[i,2]-1)*Traj[i,3]
        elseif rate_type == 11
            Wminf[i-1] = (1-Traj[i,4])*r
        elseif rate_type == 12
            Wminf[i-1] = f*Traj[i,1]*(Traj[i,1]-1)*Traj[i,4]
        else
            error("You've done something wrong here Jacob.\n")
        end
        Back -= log(Wminf[i-1])
    end
    ΔS = For + Back
    print("$For\n")
    print("$Back\n")
    print("$ΔS\n")

    plot(Tims,Traj[:,1:2])
    annotate!(0.5*Tf, 10, text("deltaS=$(ΔS)",:left))
    savefig("../Results/ExpressionLevelsHighA.png")
    plot(Traj[:,1],Traj[:,2])
    savefig("../Results/Trajs.png")
    # Don't feel that this is showing me anything particulary intresting
    #ents = log.(Wf./Wminf)
    #plot(ents)
    #savefig("../Results/EntropyProduction.png")
end

@time gillespie()
