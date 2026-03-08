# Molecular Dynamics Expert-Policy Agent

This repository is now sequence-first and LLM-policy-driven.

The objective is to optimize whether an LLM can produce the same setup decisions a human MD expert would make when given a protein sequence, then use that protocol to build a runnable OpenMM simulation workflow.

## Core Design

The repo has two complementary loops:

1. Expert-policy scoring loop (for harness optimization)
- Input: sequence task
- Policy output: expert protocol (`propose_protocol(task, pdb_candidates)`)
- Evaluator output: score for expert-step quality and simulation readiness
- No heavy MD execution required for scoring

2. Runtime execution loop (for actual setup/simulation)
- Input: sequence task + policy
- Executes sequence search, structure selection, optional sequence reconciliation, preparation, system setup
- Optional smoke test simulation
- Optional production simulation

## Key Files

- `simulation_agent/policy.py`
  - Mutatable LLM policy contract
  - Defines `propose_protocol(task, pdb_candidates)`

- `simulation_agent/setup_targets.py`
  - Sequence benchmark tasks (`TRAIN_TASKS`, `HOLDOUT_TASKS`)
  - Required expert-step tags per task

- `simulation_agent/setup_evaluator.py`
  - Harness-facing evaluator (`evaluate(program_path)`)
  - Scores schema quality, target grounding, expert-step coverage, structure strategy, setup readiness, validation checks

- `simulation_agent/protocol_executor.py`
  - Execution bridge from policy protocol to runnable setup/simulation

- `simulation_agent/agent.py`
  - Shared setup core
  - Sequence-first workflow: `run_sequence_workflow(...)`

- `simulation_agent/main.py`
  - CLI for sequence benchmark tasks and custom sequence runs

## Installation

```bash
pip install -r requirements.txt
```

## Local Policy Evaluation

Evaluate current policy against sequence expert-step benchmarks:

```bash
python -m simulation_agent.evaluate_setup_policy simulation_agent/policy.py
```

## Harness Optimization (SkyDiscover)

```bash
cd /Users/peterobi/Documents/ideas/harness/skydiscover-main
uv run skydiscover-run \
  /Users/peterobi/Documents/ideas/MD_simulation_agents/simulation_agent/policy.py \
  /Users/peterobi/Documents/ideas/MD_simulation_agents/simulation_agent/setup_evaluator.py \
  --search adaevolve \
  --model gpt-5 \
  --iterations 30 \
  --agentic
```

## Sequence Runtime Execution

Run a predefined sequence task with optional smoke test:

```bash
python -m simulation_agent.main task hba1_reference --run-smoke-test --smoke-test-steps 250
```

List available sequence tasks:

```bash
python -m simulation_agent.main task --list-tasks
```

Run a custom sequence:

```bash
python -m simulation_agent.main sequence "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF" --run-smoke-test
```

Run full production on top of setup:

```bash
python -m simulation_agent.main task adrb2_reference --run-production --production-steps 100000
```

## Notes

- Sequence reconciliation is treated as one expert decision step, not as a hardcoded WT-only rule.
- API keys are required for the external LLM optimizer loop, not for the local setup/simulation code path itself.
- Runtime execution modes require `openmm` and `pdbfixer` installed locally.
