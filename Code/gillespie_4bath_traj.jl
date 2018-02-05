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

function gillespie()
    Ti = 0.0  # initial time
    Tf = 30#000000 # end of simulation time in s
    batchsize = 10000
    Traj = zeros(0,4)
    Tims = []
    traj = zeros(batchsize,4)
    tims = zeros(batchsize)

    t = Ti # initialise
    A = convert(Int64, Ω)
    B = convert(Int64, Ω)
    W = convert(Int64, div((N - A - B), 2))
    S = convert(Int64, div((N - A - B), 2))
    if rem((N - A - B), 2) != 0
        error("Error: Inappropriate choice of initial value of A or B")
    end
    i = 1
    traj[1,:] = [ A B W S ]
    tims[1] = t
    # Main loop
    while t <= Tf
        # set prior A and Bs

        # rates
        rates = [r*k*S/(r+f*B*(B-1)), kmin*A, K*A, Kmin*W, r*q*S/(r+f*A*(A-1)), qmin*B, Q*B, Qmin*W, F]
        rs = rates/sum(rates)

        rone = rand() # first random number

        # which reaction?
        if rone < rs[1]
            A += 1
            S -= 1
        elseif rone < (rs[1]+rs[2])
            A -= 1
            S += 1
        elseif rone < (rs[1]+rs[2]+rs[3])
            A -= 1
            W += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4])
            A += 1
            W -= 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5])
            B += 1
            S -= 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6])
            B -= 1
            S += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7])
            B -= 1
            W += 1
        elseif rone < (rs[1]+rs[2]+rs[3]+rs[4]+rs[5]+rs[6]+rs[7]+rs[8])
            B += 1
            W -= 1
        else
            S += 1
            W -= 1
        end

        # update time
        rtwo = rand()  # second random number
        t = t - log(rtwo)/sum(rates)

        i += 1
        traj[i,:] = [ A B W S]
        tims[i] = t
        if i == batchsize
            Traj = vcat(Traj,traj)
            Tims = vcat(Tims,tims)
            traj = zeros(batchsize,4)
            tims = zeros(batchsize)
            i = 0
        end
    end
    # Need to now cat the final incomplete batch
    trajverytemp = fill(0, i, 4)
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
    plot(Tims,Traj)
    savefig("../Results/ExpressionLevels.png")
    # Now need to write some code that calculates the entropy production of this trajectory
    for i = 2:length(Tims)
        #calculate rate of going from pos[i] to pos[i+1]
    end
    for i = length(Tims):2
        # calculate rate of going from pos[i] to pos[i-1]
    end
end

@time gillespie()
