"""CLI entry points for sequence-first MD workflows."""

import argparse
import json

from simulation_agent.setup_targets import get_all_tasks, get_task_by_id


def _build_parser():
    """Create CLI parser for sequence-task and custom-sequence modes."""
    parser = argparse.ArgumentParser(
        description="Run sequence-first MD workflows with expert protocols."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    task_parser = subparsers.add_parser("task", help="Run a predefined sequence benchmark task.")
    task_parser.add_argument(
        "task_id",
        type=str,
        nargs="?",
        default=None,
        help="Task identifier from setup_targets.py",
    )
    task_parser.add_argument("--policy-path", type=str, default="simulation_agent/policy.py")
    task_parser.add_argument("--run-smoke-test", action="store_true")
    task_parser.add_argument("--smoke-test-steps", type=int, default=250)
    task_parser.add_argument("--run-production", action="store_true")
    task_parser.add_argument("--production-steps", type=int, default=None)
    task_parser.add_argument("--list-tasks", action="store_true")

    sequence_parser = subparsers.add_parser("sequence", help="Run a custom input sequence through the policy.")
    sequence_parser.add_argument("sequence", type=str, help="Protein sequence string")
    sequence_parser.add_argument("--task-id", type=str, default="custom_sequence")
    sequence_parser.add_argument("--task-name", type=str, default="Custom sequence task")
    sequence_parser.add_argument("--policy-path", type=str, default="simulation_agent/policy.py")
    sequence_parser.add_argument("--run-smoke-test", action="store_true")
    sequence_parser.add_argument("--smoke-test-steps", type=int, default=250)
    sequence_parser.add_argument("--run-production", action="store_true")
    sequence_parser.add_argument("--production-steps", type=int, default=None)

    return parser


def _run_task_mode(args):
    """Execute a predefined sequence task via expert policy and runtime bridge."""
    if args.list_tasks:
        print(json.dumps([task["id"] for task in get_all_tasks()], indent=2))
        return

    if not args.task_id:
        raise ValueError("task_id is required unless --list-tasks is provided.")

    from simulation_agent.protocol_executor import execute_sequence_task

    task = get_task_by_id(args.task_id)
    if task is None:
        raise ValueError(f"Unknown task_id '{args.task_id}'. Use --list-tasks to inspect valid IDs.")

    report = execute_sequence_task(
        task=task,
        policy_path=args.policy_path,
        run_smoke_test=args.run_smoke_test,
        smoke_test_steps=args.smoke_test_steps,
        run_production=args.run_production,
        production_steps=args.production_steps,
        verbose=True,
    )
    print(json.dumps(report, indent=2))


def _run_sequence_mode(args):
    """Execute a custom sequence via expert policy and runtime bridge."""
    from simulation_agent.protocol_executor import execute_sequence_task

    task = {
        "id": args.task_id,
        "name": args.task_name,
        "sequence": args.sequence,
        "metadata": {},
        "expected_target": {},
        "required_expert_steps": [
            "sequence_quality_check",
            "template_search",
            "structure_selection",
            "sequence_reconciliation",
            "protein_preparation",
            "system_parameterization",
            "solvation_and_ionization",
            "simulation_preflight_checks",
        ],
    }

    report = execute_sequence_task(
        task=task,
        policy_path=args.policy_path,
        run_smoke_test=args.run_smoke_test,
        smoke_test_steps=args.smoke_test_steps,
        run_production=args.run_production,
        production_steps=args.production_steps,
        verbose=True,
    )
    print(json.dumps(report, indent=2))


def main():
    """Dispatch CLI mode handlers."""
    parser = _build_parser()
    args = parser.parse_args()

    if args.mode == "task":
        _run_task_mode(args)
        return

    if args.mode == "sequence":
        _run_sequence_mode(args)
        return

    raise ValueError(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    main()
