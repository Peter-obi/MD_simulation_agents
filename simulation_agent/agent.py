from .tools import (
    fetch_pdb_data,
    PDBComparison,
    ProteinPreparer,
    SystemSetup,
    SimulationRunner,
)
from openmm.app import PDBFile

class SimulationAgent:
    def __init__(self, gene_name, uniprot_id, production_steps):
        self.gene_name = gene_name
        self.uniprot_id = uniprot_id
        self.production_steps = production_steps
        self.pdb_file_path = None
        self.reverted_pdb_file_path = None
        self.fixed_pdb_file_path = None
        self.system = None
        self.topology = None
        self.positions = None

    def run_workflow(self):
        """
        Runs the entire simulation workflow by calling the tools in sequence.

        The workflow consists of the following steps:
        1.  Fetch PDB data for the given gene name.
        2.  Download the PDB and FASTA files, and revert any mutations to match the
            wild-type sequence.
        3.  Prepare the protein for simulation using PDBFixer.
        4.  Set up the simulation system with a water box and ions.
        5.  Run the molecular dynamics simulation.
        """
        df = fetch_pdb_data(self.gene_name)
        print("PDB data fetched successfully.")
        print(df)

        comparison = PDBComparison(df, self.uniprot_id)
        self.pdb_file_path = comparison.download_pdb()
        comparison.download_fasta()
        self.reverted_pdb_file_path = comparison.revert_mutations()

        if self.reverted_pdb_file_path is None:
            pdb_to_prepare = self.pdb_file_path
        else:
            pdb_to_prepare = self.reverted_pdb_file_path

        preparer = ProteinPreparer(pdb_to_prepare)
        self.fixed_pdb_file_path = preparer.prepare_protein()

        pdb = PDBFile(self.fixed_pdb_file_path)
        setup = SystemSetup(pdb)
        self.system, self.topology, self.positions = setup.setup_system()

        runner = SimulationRunner(
            self.topology, self.system, self.positions, self.production_steps
        )
        runner.run_simulation()

        print("Agent workflow complete.")