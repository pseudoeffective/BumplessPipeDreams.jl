################################################################################
# BumplessPipeDreams.jl
#
# A Julia package for bumpless pipe dreams.
#
# Copyright (C) 2025 Dave Anderson, pseudoeffective.github.io
################################################################################

module BumplessPipeDreams

################################################################################
# Import
################################################################################

# Printf 
import Printf:
	@sprintf


################################################################################
# Export
################################################################################

export 
	BPD, Rothe, all_bpds, all_Kbpds, flat_bpds, top_bpds, 
	
	droop, undroop, Kdroop, unKdroop, flat_drops, top_drops, all_droops, 
	
	makeflat, isreduced, is_asm, isflat, 
	
	bpd2perm, bpd2word, bpd2asm, asm2bpd, dominant_part, 
	
	draw_bpd, print_all_bpds, print_all_Kbpds, print_flat_bpds, _draw_bpd_plots



################################################################################
# source files
################################################################################

include("bumpless_pipe_dreams.jl")
include("droops.jl")
include("asms.jl")
include("draw_bpds.jl") # must be included after `bpds.jl`

end # module BumplessPipeDreams
