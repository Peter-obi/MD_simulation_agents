"""Evaluator for expert-step MD setup planning from sequence tasks.

This module exposes evaluate(program_path), where program_path points to a
candidate policy file that defines propose_protocol(task, pdb_candidates).
"""

import importlib.util
import math
import os
import statistics
import sys
import time
import traceback

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CURRENT_DIR)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from simulation_agent.setup_targets import HOLDOUT_TASKS, TRAIN_TASKS
from simulation_agent.tools.fetch_pdb_data import fetch_pdb_data_by_sequence


LAST_REPORTS = []


def _load_policy_module(program_path):
    """Load a candidate policy module from disk."""
    spec = importlib.util.spec_from_file_location("candidate_policy", program_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import policy module from {program_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "propose_protocol"):
        raise AttributeError("Policy module must define propose_protocol(task, pdb_candidates).")
    if not callable(module.propose_protocol):
        raise TypeError("propose_protocol must be callable.")
    return module


def _safe_float(value):
    """Convert to float where possible."""
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result):
        return None
    return result


def _normalize_step_list(steps):
    """Normalize expert-step tags into a lowercase set."""
    if not isinstance(steps, list):
        return set()
    return {str(step).strip().lower() for step in steps if str(step).strip()}


def _fetch_candidates_for_task(task):
    """Fetch structure candidates for a sequence task."""
    sequence = task.get("sequence") or ""
    if not sequence:
        return None, "candidate_fetch_fail"

    try:
        df = fetch_pdb_data_by_sequence(sequence)
    except Exception:
        df = None

    if df is not None and not df.empty:
        return df, None
    return None, "candidate_fetch_fail"


def _score_schema(protocol):
    """Score protocol schema completeness."""
    if not isinstance(protocol, dict):
        return 0.0

    required_keys = {
        "target",
        "expert_steps",
        "structure_strategy",
        "system_setup",
        "simulation_protocol",
        "validation_checks",
    }
    present = sum(1 for key in required_keys if key in protocol)
    return present / len(required_keys)


def _score_target_grounding(protocol, task):
    """Score whether protocol target block correctly anchors expected identifiers."""
    target = protocol.get("target") or {}
    expected = task.get("expected_target") or {}
    if not isinstance(target, dict):
        return 0.0

    score = 0.0
    if str(target.get("gene_name", "")).upper() == str(expected.get("gene_name", "")).upper():
        score += 0.5
    if str(target.get("uniprot_id", "")).upper() == str(expected.get("uniprot_id", "")).upper():
        score += 0.5
    return score


def _score_step_coverage(protocol, task):
    """Score coverage of required expert-step tags."""
    required_steps = {
        str(step).strip().lower() for step in (task.get("required_expert_steps") or []) if str(step).strip()
    }
    if not required_steps:
        return 1.0

    provided = _normalize_step_list(protocol.get("expert_steps"))
    matched = len(required_steps.intersection(provided))
    return matched / len(required_steps)


def _score_structure_strategy(protocol, pdb_candidates):
    """Score structure-selection quality and sequence-reconciliation intent."""
    strategy = protocol.get("structure_strategy") or {}
    if not isinstance(strategy, dict):
        return 0.0

    selected_pdb_id = str(strategy.get("selected_pdb_id", "")).strip().upper()
    selected_score = 0.0
    if selected_pdb_id and pdb_candidates is not None and not pdb_candidates.empty:
        known_ids = {str(x).strip().upper() for x in pdb_candidates.get("ID", []).tolist()}
        if selected_pdb_id in known_ids:
            selected_score = 1.0

    reconcile = str(strategy.get("sequence_reconciliation", "")).strip().lower()
    reconcile_score = 1.0 if reconcile in {"required", "preferred", "optional", "reference"} else 0.0
    return 0.7 * selected_score + 0.3 * reconcile_score


def _score_setup_readiness(protocol):
    """Score whether protocol specifies runnable setup and simulation parameters."""
    system_setup = protocol.get("system_setup") or {}
    simulation_protocol = protocol.get("simulation_protocol") or {}

    if not isinstance(system_setup, dict) or not isinstance(simulation_protocol, dict):
        return 0.0

    setup_required = {
        "forcefield_files",
        "water_model",
        "padding_nm",
        "nonbonded_method",
        "constraints",
    }
    simulation_required = {"minimize", "equilibration_steps", "production_steps"}

    setup_score = sum(1 for key in setup_required if key in system_setup) / len(setup_required)
    simulation_score = (
        sum(1 for key in simulation_required if key in simulation_protocol) / len(simulation_required)
    )
    return 0.6 * setup_score + 0.4 * simulation_score


def _score_validation_quality(protocol):
    """Score preflight/validation checks proposed by the policy."""
    checks = protocol.get("validation_checks")
    if not isinstance(checks, list) or not checks:
        return 0.0

    normalized = {str(check).strip().lower() for check in checks if str(check).strip()}
    if not normalized:
        return 0.0

    desired_signals = {
        "sequence",
        "structure",
        "charge",
        "topology",
        "forcefield",
        "equilibration",
        "integrity",
    }
    overlap = 0
    for token in desired_signals:
        if any(token in check for check in normalized):
            overlap += 1

    return min(1.0, overlap / 4.0)


def _categorize_failure(component_scores, candidate_fetch_error):
    """Map low component scores to actionable failure categories."""
    if candidate_fetch_error:
        return "candidate_fetch_fail"

    if component_scores["schema_score"] < 0.5:
        return "schema_fail"
    if component_scores["target_score"] < 0.5:
        return "target_grounding_fail"
    if component_scores["step_coverage_score"] < 0.6:
        return "missing_expert_steps"
    if component_scores["setup_readiness_score"] < 0.6:
        return "setup_readiness_fail"
    return "low_quality_protocol"


def _evaluate_task(policy_module, task, split):
    """Evaluate a single sequence task and return per-task report."""
    pdb_candidates, fetch_error = _fetch_candidates_for_task(task)

    if fetch_error:
        return {
            "split": split,
            "task_id": task.get("id"),
            "task_name": task.get("name"),
            "score": 0.0,
            "success": False,
            "failure_category": fetch_error,
            "schema_score": 0.0,
            "target_score": 0.0,
            "step_coverage_score": 0.0,
            "structure_strategy_score": 0.0,
            "setup_readiness_score": 0.0,
            "validation_score": 0.0,
            "num_candidates": 0,
            "selected_pdb_id": None,
            "error": "No structure candidates available for task.",
        }

    try:
        protocol = policy_module.propose_protocol(task, pdb_candidates.copy())
    except Exception as exc:
        return {
            "split": split,
            "task_id": task.get("id"),
            "task_name": task.get("name"),
            "score": 0.0,
            "success": False,
            "failure_category": "policy_runtime_fail",
            "schema_score": 0.0,
            "target_score": 0.0,
            "step_coverage_score": 0.0,
            "structure_strategy_score": 0.0,
            "setup_readiness_score": 0.0,
            "validation_score": 0.0,
            "num_candidates": int(len(pdb_candidates)),
            "selected_pdb_id": None,
            "error": str(exc),
        }

    schema_score = _score_schema(protocol)
    target_score = _score_target_grounding(protocol, task)
    step_coverage_score = _score_step_coverage(protocol, task)
    structure_strategy_score = _score_structure_strategy(protocol, pdb_candidates)
    setup_readiness_score = _score_setup_readiness(protocol)
    validation_score = _score_validation_quality(protocol)

    combined = (
        0.15 * schema_score
        + 0.20 * target_score
        + 0.20 * step_coverage_score
        + 0.15 * structure_strategy_score
        + 0.20 * setup_readiness_score
        + 0.10 * validation_score
    )

    component_scores = {
        "schema_score": schema_score,
        "target_score": target_score,
        "step_coverage_score": step_coverage_score,
        "structure_strategy_score": structure_strategy_score,
        "setup_readiness_score": setup_readiness_score,
        "validation_score": validation_score,
    }

    failure_category = _categorize_failure(component_scores, candidate_fetch_error=None)
    success = combined >= 0.75

    structure_strategy = protocol.get("structure_strategy") if isinstance(protocol, dict) else {}
    selected_pdb_id = None
    if isinstance(structure_strategy, dict):
        selected_pdb_id = structure_strategy.get("selected_pdb_id")

    return {
        "split": split,
        "task_id": task.get("id"),
        "task_name": task.get("name"),
        "score": round(combined, 6),
        "success": bool(success),
        "failure_category": "success" if success else failure_category,
        "schema_score": round(schema_score, 6),
        "target_score": round(target_score, 6),
        "step_coverage_score": round(step_coverage_score, 6),
        "structure_strategy_score": round(structure_strategy_score, 6),
        "setup_readiness_score": round(setup_readiness_score, 6),
        "validation_score": round(validation_score, 6),
        "num_candidates": int(len(pdb_candidates)),
        "selected_pdb_id": selected_pdb_id,
        "error": None,
    }


def _aggregate_reports(reports):
    """Aggregate per-task reports into scalar metrics."""
    if not reports:
        return {
            "score": 0.0,
            "success_rate": 0.0,
            "schema_score": 0.0,
            "target_score": 0.0,
            "step_coverage_score": 0.0,
            "structure_strategy_score": 0.0,
            "setup_readiness_score": 0.0,
            "validation_score": 0.0,
            "failure_counts": {},
        }

    scores = [report["score"] for report in reports]
    schema_scores = [report["schema_score"] for report in reports]
    target_scores = [report["target_score"] for report in reports]
    step_scores = [report["step_coverage_score"] for report in reports]
    structure_scores = [report["structure_strategy_score"] for report in reports]
    readiness_scores = [report["setup_readiness_score"] for report in reports]
    validation_scores = [report["validation_score"] for report in reports]

    failure_counts = {}
    for report in reports:
        category = report.get("failure_category") or "unknown"
        if category != "success":
            failure_counts[category] = failure_counts.get(category, 0) + 1

    success_rate = sum(1 for report in reports if report.get("success")) / len(reports)

    return {
        "score": statistics.fmean(scores),
        "success_rate": success_rate,
        "schema_score": statistics.fmean(schema_scores),
        "target_score": statistics.fmean(target_scores),
        "step_coverage_score": statistics.fmean(step_scores),
        "structure_strategy_score": statistics.fmean(structure_scores),
        "setup_readiness_score": statistics.fmean(readiness_scores),
        "validation_score": statistics.fmean(validation_scores),
        "failure_counts": failure_counts,
    }


def evaluate_with_details(program_path):
    """Evaluate a policy and return aggregate metrics and per-task reports."""
    started_at = time.time()
    policy_module = _load_policy_module(program_path)

    train_reports = [_evaluate_task(policy_module, task, "train") for task in TRAIN_TASKS]
    holdout_reports = [_evaluate_task(policy_module, task, "holdout") for task in HOLDOUT_TASKS]

    train_agg = _aggregate_reports(train_reports)
    holdout_agg = _aggregate_reports(holdout_reports)

    if HOLDOUT_TASKS:
        combined_score = 0.8 * train_agg["score"] + 0.2 * holdout_agg["score"]
    else:
        combined_score = train_agg["score"]

    reports = train_reports + holdout_reports
    full_agg = _aggregate_reports(reports)

    failure_counts = full_agg["failure_counts"]
    total_tasks = len(reports)

    metrics = {
        "combined_score": round(combined_score, 6),
        "train_score": round(train_agg["score"], 6),
        "holdout_score": round(holdout_agg["score"], 6),
        "overall_success_rate": round(full_agg["success_rate"], 6),
        "train_success_rate": round(train_agg["success_rate"], 6),
        "holdout_success_rate": round(holdout_agg["success_rate"], 6),
        "avg_schema_score": round(full_agg["schema_score"], 6),
        "avg_target_score": round(full_agg["target_score"], 6),
        "avg_step_coverage_score": round(full_agg["step_coverage_score"], 6),
        "avg_structure_strategy_score": round(full_agg["structure_strategy_score"], 6),
        "avg_setup_readiness_score": round(full_agg["setup_readiness_score"], 6),
        "avg_validation_score": round(full_agg["validation_score"], 6),
        "num_tasks": float(total_tasks),
        "num_failures": float(sum(failure_counts.values())),
        "schema_failures": float(failure_counts.get("schema_fail", 0)),
        "target_grounding_failures": float(failure_counts.get("target_grounding_fail", 0)),
        "missing_expert_steps_failures": float(failure_counts.get("missing_expert_steps", 0)),
        "setup_readiness_failures": float(failure_counts.get("setup_readiness_fail", 0)),
        "candidate_fetch_failures": float(failure_counts.get("candidate_fetch_fail", 0)),
        "policy_runtime_failures": float(failure_counts.get("policy_runtime_fail", 0)),
        "evaluation_seconds": round(time.time() - started_at, 6),
    }

    if total_tasks > 0:
        metrics["schema_failure_rate"] = round(metrics["schema_failures"] / total_tasks, 6)
        metrics["target_grounding_failure_rate"] = round(
            metrics["target_grounding_failures"] / total_tasks,
            6,
        )
        metrics["missing_expert_steps_failure_rate"] = round(
            metrics["missing_expert_steps_failures"] / total_tasks,
            6,
        )
        metrics["setup_readiness_failure_rate"] = round(
            metrics["setup_readiness_failures"] / total_tasks,
            6,
        )
        metrics["candidate_fetch_failure_rate"] = round(
            metrics["candidate_fetch_failures"] / total_tasks,
            6,
        )
    else:
        metrics["schema_failure_rate"] = 0.0
        metrics["target_grounding_failure_rate"] = 0.0
        metrics["missing_expert_steps_failure_rate"] = 0.0
        metrics["setup_readiness_failure_rate"] = 0.0
        metrics["candidate_fetch_failure_rate"] = 0.0

    return metrics, reports


def evaluate(program_path):
    """SkyDiscover entry point."""
    global LAST_REPORTS

    try:
        metrics, reports = evaluate_with_details(program_path)
        LAST_REPORTS = reports
        return metrics
    except Exception:
        LAST_REPORTS = [
            {
                "split": "meta",
                "task_id": "policy_load",
                "task_name": "policy_load",
                "score": 0.0,
                "success": False,
                "failure_category": "policy_load_fail",
                "error": traceback.format_exc(),
            }
        ]
        return {
            "combined_score": 0.0,
            "train_score": 0.0,
            "holdout_score": 0.0,
            "overall_success_rate": 0.0,
            "train_success_rate": 0.0,
            "holdout_success_rate": 0.0,
            "avg_schema_score": 0.0,
            "avg_target_score": 0.0,
            "avg_step_coverage_score": 0.0,
            "avg_structure_strategy_score": 0.0,
            "avg_setup_readiness_score": 0.0,
            "avg_validation_score": 0.0,
            "num_tasks": float(len(TRAIN_TASKS) + len(HOLDOUT_TASKS)),
            "num_failures": float(len(TRAIN_TASKS) + len(HOLDOUT_TASKS)),
            "schema_failures": float(len(TRAIN_TASKS) + len(HOLDOUT_TASKS)),
            "target_grounding_failures": 0.0,
            "missing_expert_steps_failures": 0.0,
            "setup_readiness_failures": 0.0,
            "candidate_fetch_failures": 0.0,
            "policy_runtime_failures": 0.0,
            "schema_failure_rate": 1.0,
            "target_grounding_failure_rate": 0.0,
            "missing_expert_steps_failure_rate": 0.0,
            "setup_readiness_failure_rate": 0.0,
            "candidate_fetch_failure_rate": 0.0,
            "evaluation_seconds": 0.0,
        }
