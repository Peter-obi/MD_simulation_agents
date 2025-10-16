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
        """
        # 1. Fetch PDB data
        df = fetch_pdb_data(self.gene_name)
        print("PDB data fetched successfully.")
        print(df)

        # 2. Download PDB and FASTA, and revert mutations
        comparison = PDBComparison(df, self.uniprot_id)
        self.pdb_file_path = comparison.download_pdb()
        comparison.download_fasta()
        self.reverted_pdb_file_path = comparison.revert_mutations()

        if self.reverted_pdb_file_path is None:
            pdb_to_prepare = self.pdb_file_path
        else:
            pdb_to_prepare = self.reverted_pdb_file_path

        # 3. Prepare the protein
        preparer = ProteinPreparer(pdb_to_prepare)
        self.fixed_pdb_file_path = preparer.prepare_protein()

        # 4. Set up the simulation system
        pdb = PDBFile(self.fixed_pdb_file_path)
        setup = SystemSetup(pdb)
        self.system, self.topology, self.positions = setup.setup_system()

        # 5. Run the simulation
        runner = SimulationRunner(
            self.topology, self.system, self.positions, self.production_steps
        )
        runner.run_simulation()

        print("Agent workflow complete.")