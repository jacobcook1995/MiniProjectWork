#!/usr/bin/env julia
# gillespie_4bath_traj.jl
# Julia Script for Gillespie simulation of toy driven ecosystem, this now
# only calculates the histogram rather than trajectories to save computation time
# In this case 4 resivours are considered, protiens A and B, waste W, and substrate S
# This script explicitly generates trajectories so that entropy production can be found

# Add package for plotting
using Plots
import GR

# Parameters
const Ω = 30 # system size
const k = 100 # steady state for A=k/K=1
const kmin = 10.0^-20 # set all too 10.0^-20 for now
const K = k/Ω # K=k'
const Kmin = 10.0^-20
const q = 10 # steady state for B=q/Q=1
const qmin = 10.0^-20
const Q = q/Ω # Q=q'
const Qmin = 10.0^-20
const f = 10 # Promoter switching
const r = 10
const F = 250 # removal rate
const N = 2000 # number of elements in the system
const Fmin = 10.0^-30

function gillespie()
    Ti = 0.0  # initial time
    Tf = 10000#000 # end of simulation time in s
    batchsize = 100000
    Traj = zeros(0,5)
    Tims = []
    traj = zeros(batchsize,5)
    tims = zeros(batchsize)

    t = Ti # initialise
    A = convert(Int64, 76)#0)
    B = convert(Int64, 0)#750)
    W = convert(Int64, 1900)#1200)#div((N - A - B), 2))
    S = convert(Int64, 24)#50)#div((N - A - B), 2))
    trans = 0
    if rem((N - A - B), 2) != 0
        error("Error: Inappropriate choice of initial value of A or B")
    end
    i = 1
    traj[1,:] = [ A B W S trans ]
    tims[1] = t
    # Main loop
    while t <= Tf
        # set prior A and Bs

        # rates
        rates = [r*k*S/(r+f*B*(B-1)), kmin*A, K*A, Kmin*W, r*q*S/(r+f*A*(A-1)), qmin*B, Q*B, Qmin*W, F, Fmin]
        rs = rates/sum(rates)
        rone = rand() # first random number

        # which reaction?
        if rone < rs[1]
            A += 1
            S -= 1
            trans = 1
        elseif rone < (rs[1]+rs[2])
            A -= 1
            S += 1
            trans = 2
        elseif rone < (rs[1]+rs[2]+rs[3])
            A -= 1
            W += 1
            trans = 3
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4])
            A += 1
            W -= 1
            trans = 4
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5])
            B += 1
            S -= 1
            trans = 5
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6])
            B -= 1
            S += 1
            trans = 6
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7])
            B -= 1
            W += 1
            trans = 7
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8])
            B += 1
            W -= 1
            trans = 8
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8]+rs[9])
            S += 1
            W -= 1
            trans = 9
        else
            S -= 1
            W += 1
            trans = 10
        end

        # update time
        rtwo = rand()  # second random number
        t = t - log(rtwo)/sum(rates)

        i += 1
        traj[i,:] = [ A B W S trans ]
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
            Wf[i-1] = r*k*Traj[i-1,4]/(r + f*Traj[i-1,2]*(Traj[i-1,2]-1))
        elseif rate_type == 2
            Wf[i-1] = kmin*Traj[i-1,1]
        elseif rate_type == 3
            Wf[i-1] = K*Traj[i-1,1]
        elseif rate_type == 4
            Wf[i-1] = Kmin*Traj[i-1,3]
        elseif rate_type == 5
            Wf[i-1] = r*q*Traj[i-1,4]/(r + f*Traj[i-1,1]*(Traj[i-1,1]-1))
        elseif rate_type == 6
            Wf[i-1] = qmin*Traj[i-1,2]
        elseif rate_type == 7
            Wf[i-1] = Q*Traj[i-1,2]
        elseif rate_type == 8
            Wf[i-1] = Qmin*Traj[i-1,3]
        elseif rate_type == 9
            Wf[i-1] = F
        elseif rate_type == 10
            Wf[i-1] = Fmin
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
            Wminf[i-1] = r*k*Traj[i,4]/(r + f*Traj[i,2]*(Traj[i,2]-1))
        elseif rate_type == 3
            Wminf[i-1] = Kmin*Traj[i,3]
        elseif rate_type == 4
            Wminf[i-1] = K*Traj[i,1]
        elseif rate_type == 5
            Wminf[i-1] = qmin*Traj[i,2]
        elseif rate_type == 6
            Wminf[i-1] = r*q*Traj[i,4]/(r + f*Traj[i,1]*(Traj[i,1]-1))
        elseif rate_type == 7
            Wminf[i-1] = Qmin*Traj[i,3]
        elseif rate_type == 8
            Wminf[i-1] = Q*Traj[i,2]
        elseif rate_type == 9
            Wminf[i-1] = Fmin
        elseif rate_type == 10
            Wminf[i-1] = F
        else
            error("You've done something wrong here Jacob.\n")
        end
        Back -= log(Wminf[i-1])
    end
    ΔS = For + Back
    print("$For\n")
    print("$Back\n")
    print("$ΔS\n")

    plot(Tims,Traj[:,1:4])
    annotate!(0.5*Tf, 500, text("deltaS=$(ΔS)",:left))
    savefig("../Results/ExpressionLevelsHighA.png")
    # Don't feel that this is showing me anything particulary intresting
    #ents = log.(Wf./Wminf)
    #plot(ents)
    #savefig("../Results/EntropyProduction.png")
end

@time gillespie()
