"""Fetch candidate PDB metadata from RCSB using sequence queries."""

import json
from typing import Iterable

import pandas as pd
import requests


RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"


def _post_json(url, payload, timeout=30):
    """Send a POST request and return parsed JSON content."""
    response = requests.post(url, json=payload, timeout=timeout)
    response.raise_for_status()
    return response.json()


def _fetch_entry_metadata(entry_ids: Iterable[str]) -> pd.DataFrame:
    """Fetch entry-level metadata for RCSB entry identifiers."""
    entry_ids = [entry_id for entry_id in entry_ids if entry_id]
    if not entry_ids:
        return pd.DataFrame(columns=["ID", "Initial Release Date", "PubMed ID", "DOI", "Resolution"])

    query = """
    {
        entries(entry_ids: %s) {
            rcsb_id
            rcsb_accession_info {
                initial_release_date
            }
            rcsb_primary_citation {
                pdbx_database_id_PubMed
                pdbx_database_id_DOI
            }
            rcsb_entry_info {
                resolution_combined
            }
        }
    }
    """ % json.dumps(entry_ids)

    response_json = _post_json(RCSB_GRAPHQL_URL, {"query": query})
    entries = response_json.get("data", {}).get("entries") or []

    rows = []
    for entry in entries:
        resolution = None
        entry_info = entry.get("rcsb_entry_info") or {}
        resolution_combined = entry_info.get("resolution_combined")
        if isinstance(resolution_combined, list) and resolution_combined:
            resolution = resolution_combined[0]

        citation = entry.get("rcsb_primary_citation") or {}
        accession = entry.get("rcsb_accession_info") or {}

        rows.append(
            {
                "ID": entry.get("rcsb_id"),
                "Initial Release Date": accession.get("initial_release_date"),
                "PubMed ID": citation.get("pdbx_database_id_PubMed"),
                "DOI": citation.get("pdbx_database_id_DOI"),
                "Resolution": resolution,
            }
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        df["Resolution"] = pd.to_numeric(df["Resolution"], errors="coerce")
    return df


def fetch_pdb_data_by_sequence(
    sequence: str,
    max_results: int = 25,
    identity_cutoff: float = 0.9,
    evalue_cutoff: float = 1.0,
) -> pd.DataFrame:
    """Fetch PDB candidates using a protein sequence similarity query."""
    cleaned_sequence = "".join(sequence.split()).upper()
    if not cleaned_sequence:
        raise ValueError("Sequence must be non-empty.")

    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "target": "pdb_protein_sequence",
                "value": cleaned_sequence,
                "identity_cutoff": float(identity_cutoff),
                "evalue_cutoff": float(evalue_cutoff),
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": max_results},
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined",
        },
    }

    search_json = _post_json(RCSB_SEARCH_URL, query)
    polymer_entity_ids = [item.get("identifier") for item in search_json.get("result_set", [])]

    entry_ids = []
    for polymer_id in polymer_entity_ids:
        if not polymer_id:
            continue
        entry_id = str(polymer_id).split("_", 1)[0]
        if entry_id not in entry_ids:
            entry_ids.append(entry_id)

    return _fetch_entry_metadata(entry_ids)
