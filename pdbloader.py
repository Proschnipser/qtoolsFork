import os
import sys
import urllib.request
import numpy as np
import Bio
import Bio.PDB
import Bio.SeqRecord

def download_read_pdb(pdbcode, datadir, keepfile=True):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    Then it reads and returns the structure inside.
    :param pdbcode: The standard PDB ID e.g. '3ICB'
    :param datadir: The directory where the downloaded file will be saved
    :param keepfile: if False, then the downloaded file will be deleted (default: keep the downloaded file)
    :return: a Bio.PDB Structure object or None if something went wrong
    """
    pdbfilenm = download_pdb(pdbcode, datadir)
    if pdbfilenm is None:
        return None
    struct = read_pdb(pdbcode, pdbfilenm)
    if not keepfile:
        os.remove(pdbfilenm)
    return struct
        
def download_pdb(pdbcode, datadir, downloadurl="http://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
        Note that the unencrypted HTTP protocol is used by default
        to avoid spurious OpenSSL errors...
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        # all sorts of things could have gone wrong...
        print(str(err), file=sys.stderr)
        return None

def read_cif(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 
    
def read_pdb(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 

def extract_seqrecords(pdbcode, struct):
    """
    Extracts the sequence records from a Bio.PDB structure.
    :param pdbcode: the PDB ID of the structure, needed to add a sequence ID to the result
    :param struct: a Bio.PDB.Structure object
    :return: a list of Bio.SeqRecord objects
    """
    ppb = Bio.PDB.PPBuilder()
    seqrecords = []
    for i, chain in enumerate(struct.get_chains()):
        # extract and store sequences as list of SeqRecord objects
        pps = ppb.build_peptides(chain)    # polypeptides
        seq = pps[0].get_sequence() # just take the first, hope there's no chain break
        seqid = pdbcode + chain.id
        seqrec = Bio.SeqRecord.SeqRecord(seq, id=seqid, 
            description="Sequence #{}, {}".format(i+1, seqid))
        seqrecords.append(seqrec)
    return seqrecords

def extract_multi_seqrecords(pdbcodes, structures):
    """
    Extracts the sequence records from a Bio.PDB structure.
    :param pdbcode: the PDB ID of the structure, needed to add a sequence ID to the result
    :param struct: a Bio.PDB.Structure object
    :return: a list of Bio.SeqRecord objects
    """
    ppb = Bio.PDB.PPBuilder()
    seqrecords = []
    z=0
    for struct in structures:
        z+=1
        for i, chain in enumerate(struct.get_chains()):
            # extract and store sequences as list of SeqRecord objects
            pps = ppb.build_peptides(chain)    # polypeptides
            seq = pps[0].get_sequence() # just take the first, hope there's no chain break
            seqid = pdbcodes[z-1] #+"_"+ chain.id
            seqrec = Bio.SeqRecord.SeqRecord(seq, id=seqid, 
                description="Sequence #{}, {}".format(z, seqid))
            seqrecords.append(seqrec)
        
    return seqrecords

def get_calphas(struct):
    """
    Extracts the C-alpha atoms from a PDB structure.
    :param struct: A Bio.PDB.Structure object.
    :return: A list of Bio.PDB.Atom objects representing the C-alpha atoms in `struct`.
    """
    calphas = [ atom for atom in struct.get_atoms() if atom.get_fullname() == " CA " ]
    return calphas #https://stackoverflow.com/questions/10324674/parsing-a-pdb-file-in-python

def calc_distance_matrix(atoms,flatten:bool):
    """Calculate the distance matrix for a list of atoms."""
    num_atoms = len(atoms)
    
    if flatten:
        distance_mtxflat=[]
        
        for i in range(len(atoms)):
            j=i+1
            while j<len(atoms):
                if j > i:
                    distance_mtxflat.append(atoms[i] - atoms[j])# This uses the overridden '-' operator to calculate distance
                j+=1
        return distance_mtxflat
   
    distance_matrix = np.zeros((num_atoms, num_atoms))
    for i, atom1 in enumerate(atoms):
        for j, atom2 in enumerate(atoms):
            distance_matrix[i, j] = atom1 - atom2  # This uses the overridden '-' operator to calculate distance
    
        
    return distance_matrix #https://bioinformatics.stackexchange.com/questions/22675/how-can-i-compute-distance-matrices-from-a-pdb-file