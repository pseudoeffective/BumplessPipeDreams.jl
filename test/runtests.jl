using BumplessPipeDreams
using Test


if isempty(ARGS)
    tests = [ 
        "bpds.jl", 
	"draw_bpds.jl"
	]
else
    tests = ARGS
end

for test in tests
    include(test)
end

