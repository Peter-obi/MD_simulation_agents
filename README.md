# Modular Simulation Agent

This project provides a modular, agentic workflow for running molecular dynamics (MD) simulations of proteins. Given a gene name and a UniProt ID, this agent can fetch the relevant protein structure data, prepare the protein for simulation, and run a full MD simulation using OpenMM.

The codebase is designed to be modular, with each step of the workflow implemented as a separate "tool." This architecture is intended to be used by a Language Learning Model (LLM) agent to flexibly run simulation pipelines.

## Features

*   **Data Fetching:** Automatically queries the RCSB PDB database to find and download the best available crystal structure for a given gene.
*   **Sequence Comparison:** Compares the downloaded PDB sequence with the canonical sequence from UniProt and reverts any mutations to the wild type.
*   **Protein Preparation:** Uses `PDBFixer` to prepare the PDB file for simulation by adding missing atoms, hydrogens, and handling non-standard residues.
*   **System Solvation:** Sets up a simulation-ready system by solvating the protein in a water box with appropriate ions using OpenMM.
*   **MD Simulation:** Runs a standard MD simulation workflow, including energy minimization, NVT and NPT equilibration, and a final production run.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2.  **Install dependencies:**
    The required Python packages are listed in `requirements.txt`.

    First, install `pdbfixer` directly from the OpenMM GitHub repository, as it is not available on PyPI:
    ```bash
    pip install git+https://github.com/openmm/pdbfixer.git
    ```

    Then, install the rest of the dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

The main script `simulation_agent/main.py` orchestrates the entire simulation workflow. You can run it from the command line with the following arguments:

```bash
python simulation_agent/main.py <gene_name> <uniprot_id> <production_steps>
```

*   `gene_name`: The official gene symbol for the protein of interest (e.g., "OPRM1").
*   `uniprot_id`: The UniProt accession number for the protein (e.g., "P35372").
*   `production_steps`: The number of steps for the final production MD run.

### Example

To run a 10,000-step simulation for the human Mu-opioid receptor (OPRM1), use the following command:

```bash
python simulation_agent/main.py OPRM1 P35372 10000
```

The script will perform all steps from data fetching to simulation and save the output trajectory (`production.dcd`) and log file (`production.log`) in the root directory.

## Future Work

*   Integration with AlphaFold for structures without a PDB entry.
*   Advanced analysis modules using libraries like MDTraj and MDAnalysis.
*   Support for membrane protein simulations.
*   Customizable equilibration and production parameters.