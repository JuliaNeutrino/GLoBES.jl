# Simple output functions for the GLoBES examples (port of myio.c)

using Printf

const THEFILE = Ref("")

"""
    InitOutput(filename, headline)

Start a new output file and write `headline` to it. An empty filename means
standard output.
"""
function InitOutput(filename::AbstractString, headline::AbstractString)
    THEFILE[] = filename
    if isempty(filename)
        print(headline)
    else
        open(io -> print(io, headline), filename, "w")
    end
end

function output_line(line)
    if isempty(THEFILE[])
        print(line)
    else
        open(io -> print(io, line), THEFILE[], "a")
    end
end

"""
    AddToOutput(n1, n2, n3)

Add a line with three numbers to the current output.
"""
AddToOutput(n1, n2, n3) = output_line(@sprintf("%g %g %f \n", n1, n2, n3))

"""
    AddToOutput2(n1, n2)

Add a line with two numbers to the current output.
"""
AddToOutput2(n1, n2) = output_line(@sprintf("%g %g \n", n1, n2))
