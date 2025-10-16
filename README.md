# Molecular Dynamics Simulation Agent

This project provides a modular, agent-based system for running molecular dynamics (MD) simulations using OpenMM. The agent orchestrates a series of tools to fetch protein data, prepare the protein, and run a simulation, providing a flexible and extensible framework for MD simulations.

## Architecture

The project is structured around a central `SimulationAgent` that coordinates a series of tools, each responsible for a specific step in the simulation workflow:

*   **`fetch_pdb_data`**: Fetches protein data from the RCSB PDB.
*   **`PDBComparison`**: Downloads PDB and FASTA files and reverts mutations.
*   **`ProteinPreparer`**: Prepares the protein for simulation using PDBFixer.
*   **`SystemSetup`**: Sets up the simulation system with a water box and ions.
*   **`SimulationRunner`**: Runs the MD simulation.

This modular design allows for easy extension and modification of the simulation workflow.

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

The main entry point for the simulation is `simulation_agent/main.py`, which initializes and runs the `SimulationAgent`. The script takes the following command-line arguments:

*   `gene_name`: The gene name for the protein of interest (e.g., "OPRM1").
*   `uniprot_id`: The UniProt ID for the protein (e.g., "P35372").
*   `production_steps`: The number of steps for the production run.

### Example

To run a simulation for the Mu-opioid receptor (OPRM1) for 1,000,000 steps:

```bash
python -m simulation_agent.main OPRM1 P35372 1000000
```