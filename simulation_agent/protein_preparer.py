from pdbfixer import PDBFixer
from openmm.app import PDBFile

class ProteinPreparer:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.fixer = PDBFixer(filename=pdb_file)
        self.fixed_pdb_file = "fixed.pdb"

    def prepare_protein(self):
        """
        Prepares the protein for simulation by adding missing atoms and hydrogens.
        """
        self.fixer.findMissingResidues()
        self.fixer.findNonstandardResidues()
        self.fixer.replaceNonstandardResidues()
        self.fixer.removeHeterogens(False)
        self.fixer.findMissingAtoms()
        self.fixer.addMissingAtoms()
        self.fixer.addMissingHydrogens(7.0)

        with open(self.fixed_pdb_file, 'w') as f:
            PDBFile.writeFile(self.fixer.topology, self.fixer.positions, f)

        print(f"Protein prepared and saved to {self.fixed_pdb_file}")
        return self.fixed_pdb_file