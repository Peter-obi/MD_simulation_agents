from Bio.PDB import PDBList, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import requests
import warnings

def download_and_compare(df, uniprot_id):
    lowest_resolution_entry = df.loc[df['Resolution'].idxmin()]
    pdb_id = lowest_resolution_entry['ID']
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        parser = MMCIFParser()
        structure = parser.get_structure('structure', pdb_file)
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in standard_aa_names:
                    seq.append(residue.get_resname())
    pdb_seq = ''.join(seq)
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    uniprot_seq = response.text.split('\n', 1)[1].replace('\n', '')
    mutations = []
    for i, (pdb_aa, uniprot_aa) in enumerate(zip(pdb_seq, uniprot_seq), start=1):
        if pdb_aa != uniprot_aa:
            mutations.append(f"{uniprot_aa}{i}{pdb_aa}")
    print("Mutations:")
    if mutations:
        print(', '.join(mutations))
    else:
        print("No mutations found.")
    return pdb_file
