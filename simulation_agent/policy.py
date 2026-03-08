"""Baseline expert-style sequence-to-MD setup strategy policy.

SkyDiscover should mutate this file.
"""


def _rank_candidates(pdb_candidates):
    """Rank PDB candidates using resolution and release date heuristics."""
    if "Resolution" in pdb_candidates.columns:
        return pdb_candidates.sort_values(
            by=["Resolution", "Initial Release Date"],
            ascending=[True, False],
            na_position="last",
        )
    return pdb_candidates


def _build_target_block(task):
    """Build a normalized target descriptor from task metadata."""
    expected_target = task.get("expected_target") or {}
    metadata = task.get("metadata") or {}
    return {
        "gene_name": expected_target.get("gene_name"),
        "uniprot_id": expected_target.get("uniprot_id"),
        "species": metadata.get("species"),
        "family": metadata.get("family"),
    }


def propose_protocol(task, pdb_candidates):
    """Propose an expert protocol for turning a sequence into a runnable MD setup.

    Args:
        task: Benchmark task dict containing sequence and metadata.
        pdb_candidates: DataFrame from sequence-based PDB search.

    Returns:
        Protocol dictionary used by evaluator and runtime execution bridge.
    """
    if pdb_candidates is None or len(pdb_candidates) == 0:
        raise ValueError("No PDB candidates available for protocol construction.")

    ranked = _rank_candidates(pdb_candidates)
    selected_row = ranked.iloc[0]
    selected_pdb_id = str(selected_row["ID"])

    expert_steps = [
        "sequence_quality_check",
        "template_search",
        "structure_selection",
        "sequence_reconciliation",
        "protein_preparation",
        "system_parameterization",
        "solvation_and_ionization",
        "simulation_preflight_checks",
    ]

    return {
        "protocol_version": "1.0",
        "task_id": task.get("id"),
        "target": _build_target_block(task),
        "expert_steps": expert_steps,
        "structure_strategy": {
            "selected_pdb_id": selected_pdb_id,
            "selection_rule": "lowest_resolution_then_newest_release",
            "sequence_reconciliation": "preferred",
        },
        "system_setup": {
            "forcefield_files": ["amber14-all.xml", "amber14/tip3p.xml"],
            "water_model": "tip3p",
            "padding_nm": 1.0,
            "nonbonded_method": "PME",
            "nonbonded_cutoff_nm": 1.0,
            "constraints": "HBonds",
            "ewald_error_tolerance": 0.0005,
        },
        "simulation_protocol": {
            "minimize": True,
            "equilibration_steps": 1000,
            "smoke_test_steps": 250,
            "production_steps": 100000,
        },
        "validation_checks": [
            "sequence_to_structure_consistency",
            "missing_residue_review",
            "forcefield_compatibility",
            "charge_neutrality",
            "topology_integrity",
        ],
    }
