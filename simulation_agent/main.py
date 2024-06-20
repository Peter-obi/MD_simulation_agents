import argparse
from fetch_pdb_data import fetch_pdb_data
from pdb_comparison import PDBComparison

def main():
    parser = argparse.ArgumentParser(description='Fetch PDB data, download PDB file, and compare sequences.')
    parser.add_argument('gene_name', type=str, help='Gene name (e.g., "OPRM1")')
    parser.add_argument('uniprot_id', type=str, help='UniProt ID (e.g., "P35372")')
    args = parser.parse_args()
    gene_name = args.gene_name
    uniprot_id = args.uniprot_id

    df = fetch_pdb_data(gene_name)
    print(df)

    comparison = PDBComparison(df, uniprot_id)
    comparison.download_pdb()
    comparison.download_fasta()
    comparison.revert_mutations()

if __name__ == '__main__':
    main()