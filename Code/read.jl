#!/usr/bin/env julia
# read.jl
# A script to read in my file and find distinct points
using Plots
import GR

function main()
    # Assign the first command line argument to a variable called input_file
    input_file1 = "../Results/$(ARGS[1]).csv"
    input_file2 = "../Results/$(ARGS[2]).csv"
    points1 = Array{Float64}(0,2)
    points2 = Array{Float64}(0,2)
    # Open the input file for reading and close automatically at end
    open(input_file1, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for i = 1:L
                if line[i] == ','
                    comma = i
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points1 = vcat(points1, [ A B ])
        end
    end
    # now do second file
    open(input_file2, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)
            # parse line by finding commas
            comma = 0
            L = length(line)
            for i = 1:L
                if line[i] == ','
                    comma = i
                end
            end
            A = parse(Float64, line[1:(comma - 1)])
            B = parse(Float64, line[(comma + 1):L])
            points2 = vcat(points2, [ A B ])
        end
    end
    p = plot(points1[:,1],points1[:,2],label="MAP")
    p = plot!(p,points2[:,1],points2[:,2],xaxis = "A",yaxis = "B",label="gMAP")
    p = scatter!(p, [points1[1,1]], [points1[1,2]], seriescolor = :green,label="Start")
    p = scatter!(p, [points1[end,1]], [points1[end,2]], seriescolor = :red,label="Finish")
    savefig("../Results/CombinedGraph.png")
end

@time main()
