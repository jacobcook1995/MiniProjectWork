#!/usr/bin/env julia
# read3bath.jl
# A script to read in my file and find distinct points
using Plots
import GR

function main()
    # create array to hold read in data
    points1 = Array{Float64}(0,3)
    points2 = Array{Float64}(0,3)
    # check if file is provided then read this data
    if length(ARGS) > 1
        println("Reading in $(ARGS[1])")
        open(ARGS[1], "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                comma = fill(0,2)
                j = 0
                L = length(line)
                for i = 1:L
                    if line[i] == ','
                        j += 1
                        comma[j] = i
                    end
                end
                A = parse(Float64, line[1:(comma[1] - 1)])
                B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
                S = parse(Float64, line[(comma[2] + 1):L])
                points1 = vcat(points1, [ A B S ])
            end
        end
        println("Reading in $(ARGS[2])")
        open(ARGS[2], "r") do in_file
            # Use a for loop to process the rows in the input file one-by-one
            for line in eachline(in_file)
                # parse line by finding commas
                comma = fill(0,2)
                j = 0
                L = length(line)
                for i = 1:L
                    if line[i] == ','
                        j += 1
                        comma[j] = i
                    end
                end
                A = parse(Float64, line[1:(comma[1] - 1)])
                B = parse(Float64, line[(comma[1] + 1):(comma[2] - 1)])
                S = parse(Float64, line[(comma[2] + 1):L])
                points2 = vcat(points2, [ A B S ])
            end
        end
    else
        println("Need to provide 2 files to read")
        quit()
    end
    # find NM for each data set
    NG1 = size(points1,1)
    NG2 = size(points2,1)
    if NG1 != NG2
        println("Files unequal length")
        quit()
    end
    # think firstly I'd like to try and plot the path in 3D
    plot3d(points1[:,1], points1[:,2], points1[:,3], xlab = "A", ylab = "B", zlab = "S", title = "Path 1")
    savefig("../Results/Graph3D1.png")
    plot3d(points2[:,1], points2[:,2], points2[:,3], xlab = "A", ylab = "B", zlab = "S", title = "Path 2")
    savefig("../Results/Graph3D2.png")
    # make comnbined vector so both can be plotted
    points = zeros(NG1,3,2)
    for j = 1:3
        for i = 1:NG1
            points[i,j,1] = points1[i,j]
            points[i,j,2] = points2[i,j]
        end
    end
    plot3d(points[:,1,:], points[:,2,:], points[:,3,:], xlab = "A", ylab = "B", zlab = "S", title = "Path 1&2")
    savefig("../Results/Graph3D.png")
end

@time main()
