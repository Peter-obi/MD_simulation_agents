"""Simulation-agent orchestration for sequence-first MD setup and execution."""

import time
from typing import Callable

from openmm.app import PDBFile

from .tools.fetch_pdb_data import fetch_pdb_data_by_sequence
from .tools.pdb_comparison import PDBComparison
from .tools.protein_preparer import ProteinPreparer
from .tools.simulation_runner import SimulationRunner
from .tools.system_setup import SystemSetup


class SimulationAgent:
    """Coordinate structure selection, setup, and optional simulation runs from sequence."""

    def __init__(self, production_steps=100000):
        """Initialize agent runtime state."""
        self.production_steps = production_steps
        self.pdb_file_path = None
        self.reconciled_pdb_file_path = None
        self.fixed_pdb_file_path = None
        self.system = None
        self.topology = None
        self.positions = None

    def _extract_decision_fields(self, decision):
        """Normalize policy output into setup-specific fields."""
        decision = decision or {}
        if not isinstance(decision, dict):
            raise TypeError("Policy function must return a dictionary.")

        structure_strategy = decision.get("structure_strategy") or {}
        if not isinstance(structure_strategy, dict):
            raise TypeError("'structure_strategy' must be a dictionary.")

        selected_pdb_id = decision.get("selected_pdb_id")
        if selected_pdb_id is None:
            selected_pdb_id = structure_strategy.get("selected_pdb_id")

        system_setup_options = decision.get("system_setup") or {}
        if not isinstance(system_setup_options, dict):
            raise TypeError("'system_setup' must be a dictionary when provided.")

        reconcile_mode = structure_strategy.get("sequence_reconciliation", "optional")
        if not isinstance(reconcile_mode, str):
            reconcile_mode = "optional"

        expert_steps = decision.get("expert_steps") or []
        if not isinstance(expert_steps, list):
            expert_steps = []

        simulation_protocol = decision.get("simulation_protocol") or {}
        if not isinstance(simulation_protocol, dict):
            simulation_protocol = {}

        return {
            "selected_pdb_id": selected_pdb_id,
            "system_setup": system_setup_options,
            "sequence_reconciliation": reconcile_mode,
            "expert_steps": expert_steps,
            "simulation_protocol": simulation_protocol,
        }

    def _run_setup_core(self, df, reference_sequence, decision, verbose):
        """Execute the shared PDB->prepare->system setup core."""
        started_at = time.time()
        failure_stage = "unknown"

        normalized = self._extract_decision_fields(decision)
        report = {
            "success": False,
            "mode": "sequence",
            "selected_pdb_id": None,
            "selected_resolution": None,
            "num_atoms": 0,
            "num_residues": 0,
            "num_chains": 0,
            "num_particles": 0,
            "elapsed_seconds": 0.0,
            "failure_stage": None,
            "error": None,
            "expert_steps": normalized["expert_steps"],
            "sequence_reconciliation": normalized["sequence_reconciliation"],
            "simulation_protocol": normalized["simulation_protocol"],
        }

        try:
            failure_stage = "pdb_download"
            comparison = PDBComparison(df, reference_sequence=reference_sequence)
            self.pdb_file_path = comparison.download_pdb(pdb_id=normalized["selected_pdb_id"])

            selected_entry = comparison.selected_entry or {}
            report["selected_pdb_id"] = selected_entry.get("ID")
            report["selected_resolution"] = selected_entry.get("Resolution")

            failure_stage = "sequence_reconciliation"
            reconcile_mode = normalized["sequence_reconciliation"].lower()
            if reconcile_mode in {"required", "preferred", "optional", "reference"}:
                reconciled_file = comparison.reconcile_to_reference(reference_sequence=reference_sequence)
                self.reconciled_pdb_file_path = reconciled_file
            else:
                self.reconciled_pdb_file_path = None

            pdb_to_prepare = self.reconciled_pdb_file_path or self.pdb_file_path

            failure_stage = "prepare"
            preparer = ProteinPreparer(pdb_to_prepare)
            self.fixed_pdb_file_path = preparer.prepare_protein()

            failure_stage = "system_setup"
            pdb = PDBFile(self.fixed_pdb_file_path)
            setup = SystemSetup(pdb, options=normalized["system_setup"])
            self.system, self.topology, self.positions = setup.setup_system()

            report["num_atoms"] = self.topology.getNumAtoms()
            report["num_residues"] = self.topology.getNumResidues()
            report["num_chains"] = self.topology.getNumChains()
            report["num_particles"] = self.system.getNumParticles()
            report["success"] = True
            report["failure_stage"] = None
        except Exception as exc:
            report["error"] = str(exc)
            report["failure_stage"] = failure_stage
        finally:
            report["elapsed_seconds"] = round(time.time() - started_at, 3)
            if verbose:
                print(report)

        return report

    def run_sequence_workflow(
        self,
        sequence: str,
        policy_fn: Callable | None = None,
        run_smoke_test: bool = False,
        smoke_test_steps: int = 250,
        verbose: bool = True,
        task_metadata: dict | None = None,
    ):
        """Run sequence-first setup and optional smoke simulation execution."""
        cleaned_sequence = "".join(sequence.split()).upper()
        if not cleaned_sequence:
            raise ValueError("Sequence must be non-empty.")

        df = fetch_pdb_data_by_sequence(cleaned_sequence)
        if df.empty:
            return {
                "success": False,
                "mode": "sequence",
                "failure_stage": "fetch",
                "error": "No sequence-matched PDB entries were found.",
                "elapsed_seconds": 0.0,
                "smoke_test_ran": False,
                "smoke_test_success": False,
            }

        decision = {}
        if policy_fn is not None:
            decision = policy_fn(
                {
                    "sequence": cleaned_sequence,
                    "task_metadata": task_metadata or {},
                },
                df.copy(),
            )

        report = self._run_setup_core(
            df=df,
            reference_sequence=cleaned_sequence,
            decision=decision,
            verbose=verbose,
        )

        report["smoke_test_ran"] = False
        report["smoke_test_success"] = False
        report["smoke_test_steps"] = 0

        if run_smoke_test and report.get("success"):
            smoke_steps = int(smoke_test_steps)
            protocol_steps = report.get("simulation_protocol", {}).get("smoke_test_steps")
            if protocol_steps is not None:
                smoke_steps = int(protocol_steps)

            if smoke_steps > 0:
                report["smoke_test_ran"] = True
                report["smoke_test_steps"] = smoke_steps
                try:
                    runner = SimulationRunner(self.topology, self.system, self.positions, smoke_steps)
                    runner.run_simulation()
                    report["smoke_test_success"] = True
                except Exception as exc:
                    report["smoke_test_success"] = False
                    report["smoke_test_error"] = str(exc)

        return report
