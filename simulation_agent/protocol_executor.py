"""Execution bridge from sequence tasks to runnable MD setup/simulation workflows."""

import importlib.util

from simulation_agent.agent import SimulationAgent
from simulation_agent.tools.simulation_runner import SimulationRunner


def load_policy_module(program_path):
    """Load a policy module that defines propose_protocol(task, pdb_candidates)."""
    spec = importlib.util.spec_from_file_location("policy_module", program_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load policy module from {program_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "propose_protocol"):
        raise AttributeError("Policy module must define propose_protocol(task, pdb_candidates).")
    if not callable(module.propose_protocol):
        raise TypeError("propose_protocol must be callable.")
    return module


def execute_sequence_task(
    task,
    policy_path="simulation_agent/policy.py",
    run_smoke_test=True,
    smoke_test_steps=250,
    run_production=False,
    production_steps=None,
    verbose=True,
):
    """Execute setup (and optional simulation) for a benchmark task using a policy."""
    if not isinstance(task, dict):
        raise TypeError("task must be a dictionary.")
    if not task.get("sequence"):
        raise ValueError("task must include a non-empty 'sequence'.")

    policy_module = load_policy_module(policy_path)

    def policy_fn(target, pdb_candidates):
        task_payload = dict(task)
        task_payload["sequence"] = target.get("sequence")
        if target.get("task_metadata") is not None:
            task_payload["metadata"] = target.get("task_metadata")
        return policy_module.propose_protocol(task_payload, pdb_candidates)

    agent = SimulationAgent(production_steps=0)
    report = agent.run_sequence_workflow(
        sequence=task["sequence"],
        policy_fn=policy_fn,
        run_smoke_test=run_smoke_test,
        smoke_test_steps=smoke_test_steps,
        verbose=verbose,
        task_metadata=task.get("metadata"),
    )

    report["production_ran"] = False
    report["production_success"] = False
    report["production_steps"] = 0

    if run_production and report.get("success"):
        steps = production_steps
        if steps is None:
            protocol_steps = report.get("simulation_protocol", {}).get("production_steps")
            steps = int(protocol_steps) if protocol_steps is not None else 100000

        if steps > 0:
            runner = SimulationRunner(agent.topology, agent.system, agent.positions, steps)
            try:
                runner.run_simulation()
                report["production_ran"] = True
                report["production_success"] = True
                report["production_steps"] = int(steps)
            except Exception as exc:
                report["production_ran"] = True
                report["production_success"] = False
                report["production_steps"] = int(steps)
                report["production_error"] = str(exc)

    return report
