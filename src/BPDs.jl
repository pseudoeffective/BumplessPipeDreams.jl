################################################################################
# BPDs.jl
#
# A Julia package for bumpless pipe dreams.
#
# Copyright (C) 2025 Dave Anderson, pseudoeffective.github.io
################################################################################

module BPDs

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
	BPD, Rothe, all_bpds, all_Kbpds, flat_bpds, droop, undroop, Kdroop, unKdroop, 
	
	flat_drops, all_droops, makeflat, isreduced, bpd2perm, bpd2word, dominant_part, 
	
	bpd2asm, asm2bpd, is_asm, draw_bpd, print_all_bpds, print_all_Kbpds, print_flat_bpds, 
	
	_draw_bpd_plots



################################################################################
# source files
################################################################################

include("bumpless_pipe_dreams.jl")
include("asms.jl")
include("draw_bpds.jl") # must be included after `bpds.jl`

end # module BPDs
