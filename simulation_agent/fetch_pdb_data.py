import requests
import json
import pandas as pd

def fetch_pdb_data(gene_name):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "label": "full_text",
            "service": "full_text",
            "parameters": {
                "value": gene_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 25
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }
    response = requests.post(url, json=query)
    data = response.json()
    pdb_ids = [item["identifier"] for item in data.get("result_set", [])]
    url = 'https://data.rcsb.org/graphql'
    query = '''
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
    ''' % json.dumps(pdb_ids)
    response = requests.post(url, json={'query': query})
    data = response.json()
    entries = data['data']['entries']
    entry_data = []
    for entry in entries:
        entry_data.append({
            'ID': entry['rcsb_id'],
            'Initial Release Date': entry['rcsb_accession_info']['initial_release_date'],
            'PubMed ID': entry['rcsb_primary_citation']['pdbx_database_id_PubMed'],
            'DOI': entry['rcsb_primary_citation']['pdbx_database_id_DOI'],
            'Resolution': entry['rcsb_entry_info']['resolution_combined']
        })
    df = pd.DataFrame(entry_data)
    return df
