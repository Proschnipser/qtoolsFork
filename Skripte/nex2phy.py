from Bio import AlignIO
import sys

input=sys.argv[1]
output=input.replace(".nex",".phy")
input_handle = open(input, "r")
output_handle = open(output, "w")

# Format can be 'phylip' (interleaved) or 'phylip-sequential'
AlignIO.convert(input_handle, "nexus", output_handle, "phylip-relaxed")

input_handle.close()
output_handle.close()