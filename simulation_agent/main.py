import argparse
from .agent import SimulationAgent

def main():
    parser = argparse.ArgumentParser(description='Run a protein simulation workflow.')
    parser.add_argument('gene_name', type=str, help='Gene name (e.g., "OPRM1")')
    parser.add_argument('uniprot_id', type=str, help='UniProt ID (e.g., "P35372")')
    parser.add_argument('production_steps', type=int, help='Number of steps for the production run')
    args = parser.parse_args()

    agent = SimulationAgent(
        gene_name=args.gene_name,
        uniprot_id=args.uniprot_id,
        production_steps=args.production_steps,
    )
    agent.run_workflow()

if __name__ == '__main__':
    main()