# Molecular Dynamics Simulation Agent

This project is a modular, agentic system that uses a Large Language Model (LLM) to run molecular dynamics (MD) simulations using OpenMM. Given a UniProt ID and a gene name, the agent can fetch protein structures, prepare them for simulation, and run the simulation for a specified number of steps.

## Features

*   Fetches protein data from the PDB based on a gene name.
*   Downloads PDB and FASTA files.
*   Reverts mutations in the PDB file to match the wild-type sequence from UniProt.
*   Prepares the protein for simulation using PDBFixer (adds missing atoms, hydrogens, and solvent).
*   Sets up and runs a simulation using OpenMM.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2.  **Create and activate a virtual environment (recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **Install the required packages:**
    The dependencies are listed in `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```
    *Note: `pdbfixer` is installed directly from its GitHub repository to ensure the latest version is used.*

## Usage

The main script for running the simulation workflow is `simulation_agent/main.py`. It takes the following command-line arguments:

*   `gene_name`: The gene name for the protein of interest (e.g., "OPRM1").
*   `uniprot_id`: The UniProt ID for the protein (e.g., "P35372").
*   `production_steps`: The number of steps for the production run.

### Example

To run a simulation for the Mu-opioid receptor (OPRM1) for 1,000,000 steps:

```bash
python simulation_agent/main.py OPRM1 P35372 1000000
```