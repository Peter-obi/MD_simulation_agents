"""CLI utility to evaluate expert-sequence MD setup policy quality."""

import argparse
import json

from simulation_agent.setup_evaluator import evaluate_with_details


def main():
    """Run expert-policy evaluation from the command line."""
    parser = argparse.ArgumentParser(
        description="Evaluate an expert MD setup policy from sequence benchmark tasks."
    )
    parser.add_argument(
        "program_path",
        nargs="?",
        default="simulation_agent/policy.py",
        help="Path to a policy file that defines propose_protocol(task, pdb_candidates).",
    )
    args = parser.parse_args()

    metrics, reports = evaluate_with_details(args.program_path)
    print(json.dumps({"metrics": metrics, "reports": reports}, indent=2))


if __name__ == "__main__":
    main()
