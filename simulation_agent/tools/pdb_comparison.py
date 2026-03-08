"""PDB selection and sequence-reconciliation utilities."""

from pathlib import Path

from Bio import PDB
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.Data.PDBData import protein_letters_1to3
from Bio.PDB import PDBIO, PDBList


class PDBComparison:
    """Handle PDB selection and sequence-aware residue reconciliation."""

    def __init__(self, df, reference_sequence=None):
        """Initialize with candidate entries and optional reference sequence."""
        self.df = df
        self.reference_sequence = reference_sequence
        self.pdb_file = None
        self.selected_entry = None

    def _select_entry(self, pdb_id=None):
        """Select either an explicit PDB ID or the best-resolution candidate."""
        if self.df.empty:
            raise ValueError("PDB candidate table is empty.")

        if pdb_id is not None:
            matches = self.df[self.df["ID"] == pdb_id]
            if matches.empty:
                raise ValueError(f"PDB ID '{pdb_id}' was not found in the candidate table.")
            return matches.iloc[0]

        if self.df["Resolution"].notna().any():
            return self.df.loc[self.df["Resolution"].idxmin()]

        return self.df.iloc[0]

    def download_pdb(self, pdb_id=None):
        """Download the selected PDB file and store its metadata."""
        selected_entry = self._select_entry(pdb_id=pdb_id)
        self.selected_entry = selected_entry.to_dict()
        pdb_id = selected_entry["ID"]
        pdbl = PDBList()
        self.pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
        return self.pdb_file

    def _extract_pdb_sequence(self):
        """Return PDB sequence and residue handles in alignment order."""
        if not self.pdb_file:
            raise ValueError("No PDB file available. Call download_pdb() first.")

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.pdb_file)

        ppb = PDB.PPBuilder()
        pdb_seq = ""
        pdb_residues = []
        for peptide in ppb.build_peptides(structure):
            pdb_seq += str(peptide.get_sequence())
            for residue in peptide:
                pdb_residues.append(residue)

        return pdb_seq, pdb_residues

    def get_reference_sequence(self, reference_sequence=None):
        """Resolve the reference sequence from explicit input."""
        explicit = reference_sequence or self.reference_sequence
        if explicit:
            return "".join(str(explicit).split()).upper()
        return None

    def find_mismatches(self, reference_sequence=None):
        """Find residue mismatches between structure sequence and reference sequence."""
        target_sequence = self.get_reference_sequence(reference_sequence=reference_sequence)
        if not target_sequence:
            raise ValueError("A reference sequence is required to reconcile residues.")

        pdb_seq, pdb_residues = self._extract_pdb_sequence()
        alignments = pairwise2.align.localds(
            pdb_seq,
            target_sequence,
            substitution_matrices.load("BLOSUM62"),
            -12,
            -1,
        )
        if not alignments:
            return []

        aligned_pdb_seq, aligned_ref_seq = alignments[0][:2]

        mismatches = []
        pdb_index = 0
        for i in range(len(aligned_pdb_seq)):
            pdb_res = aligned_pdb_seq[i]
            ref_res = aligned_ref_seq[i]
            if pdb_res != "-":
                pdb_index += 1
            if pdb_res != "-" and ref_res != "-" and pdb_res != ref_res:
                mismatches.append((pdb_residues[pdb_index - 1], pdb_res, ref_res))

        return mismatches

    def reconcile_to_reference(self, reference_sequence=None, output_suffix="_reconciled"):
        """Rewrite mismatched residues to match the selected reference sequence."""
        mismatches = self.find_mismatches(reference_sequence=reference_sequence)
        if not mismatches:
            return None

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", self.pdb_file)

        reconciled = 0
        for residue, _pdb_aa, ref_aa in mismatches:
            ref_aa = ref_aa.upper()
            if ref_aa in protein_letters_1to3:
                residue.resname = protein_letters_1to3[ref_aa]
                reconciled += 1

        if reconciled == 0:
            return None

        stem = Path(self.pdb_file).stem
        reconciled_file = f"{stem}{output_suffix}.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(reconciled_file)
        return reconciled_file
