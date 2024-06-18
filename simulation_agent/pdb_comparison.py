from Bio import PDB
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB import PDBList, PDBParser, PDBIO
from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import BiopythonDeprecationWarning
import requests
import warnings
import traceback

class PDBComparison:
    def __init__(self, df, uniprot_id):
        self.df = df
        self.uniprot_id = uniprot_id
        self.pdb_file = None

    def download_pdb(self):
        if self.df['Resolution'].notna().any():
            lowest_resolution_entry = self.df.loc[self.df['Resolution'].idxmin()]
        else:
            print("No valid resolution values found. Selecting the first entry.")
            lowest_resolution_entry = self.df.iloc[0]

        pdb_id = lowest_resolution_entry['ID']
        pdbl = PDBList()
        self.pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        print(f"PDB file {pdb_id} downloaded successfully.")

    def download_fasta(self):
        url = f"https://www.uniprot.org/uniprot/{self.uniprot_id}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"{self.uniprot_id}.fasta", "w") as f:
                f.write(response.text)
            print(f"FASTA file {self.uniprot_id} downloaded successfully.")
        else:
            print(f"Failed to download FASTA file {self.uniprot_id}.")

    def find_mutations(self):
        # Parse the PDB file
        parser = PDB.PDBParser()
        structure = parser.get_structure("structure", self.pdb_file)

        # Extract the sequence from the PDB file
        ppb = PDB.PPBuilder()
        pdb_seq = ""
        pdb_residues = []
        for pp in ppb.build_peptides(structure):
            pdb_seq += str(pp.get_sequence())
            for res in pp:
                pdb_residues.append(res)

        # Read the FASTA file
        fasta_seq = str(SeqIO.read(f"{self.uniprot_id}.fasta", "fasta").seq)

        # Perform local alignment using the Waterman-Smith-Beyer algorithm with BLOSUM62
        alignments = pairwise2.align.localds(pdb_seq, fasta_seq, substitution_matrices.load("BLOSUM62"), -12, -1)
        aligned_pdb_seq, aligned_fasta_seq = alignments[0][:2]

        # Find mutations
        mutations = []
        pdb_index = fasta_index = 0
        for i in range(len(aligned_pdb_seq)):
            pdb_res = aligned_pdb_seq[i]
            fasta_res = aligned_fasta_seq[i]
            if pdb_res != "-":
                pdb_index += 1
            if fasta_res != "-":
                fasta_index += 1
            if pdb_res != "-" and fasta_res != "-" and pdb_res != fasta_res:
                mutations.append((pdb_residues[pdb_index - 1], pdb_res, fasta_res))

        return mutations

    def revert_mutations(self):
        mutations = self.find_mutations()
        if not mutations:
            print("No mutations found.")
            return

        print("Mutations found:")
        for mutation in mutations:
            print(f"Residue: {mutation[0].get_id()[1]}, PDB AA: {mutation[1]}, FASTA AA: {mutation[2]}")

        print("Reverting mutations...")
        parser = PDB.PDBParser()
        structure = parser.get_structure("structure", self.pdb_file)

        for mutation in mutations:
            residue = mutation[0]
            fasta_aa = mutation[2]
            residue.resname = one_to_three(fasta_aa)

        pdb_id = structure.get_id()
        io = PDBIO()
        io.set_structure(structure)
        io.save(f"{pdb_id}_reverted.pdb")
        print(f"Reverted structure saved as {pdb_id}_reverted.pdb")