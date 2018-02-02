#!/usr/bin/env julia
# read.jl
# A script to read in my file and find distinct points
function main()
    # Assign the first command line argument to a variable called input_file
    input_file = ARGS[1]
    # Assign the second command line argument to a variable called output_file
    output_file = ARGS[2]

    # Open the output file for writing
    out_file = open(output_file, "w")
    # Open the input file for reading and close automatically at end
    open(input_file, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        points = Array{Float64}(0,4) # empty array for the points to be input to
        n = Array{Int64}(0) # empty vector to store numbers
        for line in eachline(in_file)
            # Write the row of data to the output file
            # Check if it is a real point
            real = true
            distinct = true
            condition1 = false
            condition2 = false
            condition3 = false
            condition4 = false
            if line[1] == 'N'
                real = false
            end
            # if so then continue
            if real == true
                # parse line by finding commas
                comma = fill(0, 3)
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
                W = parse(Float64, line[(comma[2] + 1):(comma[3] - 1)])
                S = parse(Float64, line[(comma[3] + 1):L])

                for i = 1:size(points, 1)
                    # Now need to check if it is an orginal point
                    condition1 = A < (points[i,1] - 10.0^-5) || A > (points[i,1] + 10.0^-5)
                    condition2 = B < (points[i,2] - 10.0^-5) || B > (points[i,2] + 10.0^-5)
                    condition3 = W < (points[i,3] - 10.0^-5) || W > (points[i,3] + 10.0^-5)
                    condition4 = S < (points[i,4] - 10.0^-5) || S > (points[i,4] + 10.0^-5)
                    if condition1 != true && condition2 != true && condition3 != true && condition4 != true
                        distinct = false
                        points[i,1] = (n[i]*points[i,1] + A)/(n[i] + 1)
                        points[i,2] = (n[i]*points[i,2] + B)/(n[i] + 1)
                        points[i,3] = (n[i]*points[i,3] + W)/(n[i] + 1)
                        points[i,4] = (n[i]*points[i,4] + S)/(n[i] + 1)
                        n[i] += 1
                    end
                end
                if distinct == true
                    # if so add point to array that orginality is checked relative too
                    points = vcat(points, [ A B W S ])
                    n = vcat(n, 1)
                end
            end
        end
        for i = 1:size(points,1)
            write(out_file, "$(points[i,1]),$(points[i,2]),$(points[i,3]),$(points[i,4])\n")
        end
    end
    # Close the output file handle
    close(out_file)
end

@time main()
