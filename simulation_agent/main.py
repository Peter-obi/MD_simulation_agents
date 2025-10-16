import argparse
from fetch_pdb_data import fetch_pdb_data
from pdb_comparison import PDBComparison
from protein_preparer import ProteinPreparer
from system_setup import SystemSetup
from simulation_runner import SimulationRunner
from openmm.app import PDBFile

def main():
    parser = argparse.ArgumentParser(description='Run a protein simulation workflow.')
    parser.add_argument('gene_name', type=str, help='Gene name (e.g., "OPRM1")')
    parser.add_argument('uniprot_id', type=str, help='UniProt ID (e.g., "P35372")')
    parser.add_argument('production_steps', type=int, help='Number of steps for the production run')
    args = parser.parse_args()

    # 1. Fetch PDB data
    df = fetch_pdb_data(args.gene_name)
    print("PDB data fetched successfully.")
    print(df)

    # 2. Download PDB and FASTA, and revert mutations
    comparison = PDBComparison(df, args.uniprot_id)
    pdb_file_path = comparison.download_pdb()
    comparison.download_fasta()
    reverted_pdb_file_path = comparison.revert_mutations()

    if reverted_pdb_file_path is None:
        pdb_to_prepare = pdb_file_path
    else:
        pdb_to_prepare = reverted_pdb_file_path

    # 3. Prepare the protein
    preparer = ProteinPreparer(pdb_to_prepare)
    fixed_pdb_file_path = preparer.prepare_protein()

    # 4. Set up the simulation system
    pdb = PDBFile(fixed_pdb_file_path)
    setup = SystemSetup(pdb)
    system, topology, positions = setup.setup_system()

    # 5. Run the simulation
    runner = SimulationRunner(topology, system, positions, args.production_steps)
    runner.run_simulation()

if __name__ == '__main__':
    main()