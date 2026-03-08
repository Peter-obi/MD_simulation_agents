"""Sequence-centric benchmark tasks for expert-step MD setup planning."""

TRAIN_TASKS = [
    {
        "id": "hba1_reference",
        "name": "Human hemoglobin alpha",
        "sequence": "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "metadata": {
            "species": "Homo sapiens",
            "family": "globin",
            "notes": "Stable soluble protein commonly used for structural benchmarking.",
        },
        "expected_target": {
            "gene_name": "HBA1",
            "uniprot_id": "P69905",
        },
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
    },
    {
        "id": "adrb2_reference",
        "name": "Beta-2 adrenergic receptor",
        "sequence": "MGQPGNGSAFLLAPNGSHAPDHDVTQQRDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFW",
        "metadata": {
            "species": "Homo sapiens",
            "family": "GPCR",
            "notes": "Membrane receptor; expert setup should account for membrane-relevant choices.",
        },
        "expected_target": {
            "gene_name": "ADRB2",
            "uniprot_id": "P07550",
        },
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
    },
    {
        "id": "egfr_variant_like",
        "name": "EGFR variant-like sequence",
        "sequence": "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTNNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIAHNQVRQVPLQRSL",
        "metadata": {
            "species": "Homo sapiens",
            "family": "receptor tyrosine kinase",
            "notes": "Expert workflow should explicitly reason about sequence reconciliation and mutation handling options.",
        },
        "expected_target": {
            "gene_name": "EGFR",
            "uniprot_id": "P00533",
        },
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
    },
]

HOLDOUT_TASKS = [
    {
        "id": "braf_variant_like",
        "name": "BRAF variant-like sequence",
        "sequence": "MAALSGGGGGGAEPGQALFNGDMEPEAGAGPGSLTQPEAAPRLPQHPRAAASLAPNFGQPRQSVLRGQQAATNADQENKRQGQGQTSLQESDGSQGQPGQGQGS",
        "metadata": {
            "species": "Homo sapiens",
            "family": "kinase",
            "notes": "Holdout task to test generalization of expert-step planning.",
        },
        "expected_target": {
            "gene_name": "BRAF",
            "uniprot_id": "P15056",
        },
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
    },
]


def get_all_tasks():
    """Return all benchmark tasks."""
    return TRAIN_TASKS + HOLDOUT_TASKS


def get_task_by_id(task_id: str):
    """Return a task by identifier or None when it is not defined."""
    for task in get_all_tasks():
        if task.get("id") == task_id:
            return task
    return None
