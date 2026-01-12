from pdbloader import *
import os

#struct= download_read_pdb("7R3M", os.getcwd()+"")
struct= read_pdb("7R3M","7R3M.pdb")
calphas= get_calphas(struct)
print(extract_seqrecords("7R3M",struct))
print(struct)


atoms = get_calphas(struct)
print(calc_distance_matrix(atoms))