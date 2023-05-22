import requests
import json
import pandas as pd

# Change the gene name to the desired gene_name (this can be gotten from uniprot)
gene_name = "OPRM1"

url = " https://search.rcsb.org/rcsbsearch/v2/query"

# Define the query
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


# Make the request
response = requests.post(url, json=query)

# Extract the data from the response
data = response.json()

# Get the PDB_ids
pdb_ids = [item["identifier"] for item in data.get("result_set", [])]

url = 'https://data.rcsb.org/graphql'

# Define the query
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

# Make the request
response = requests.post(url, json={'query': query})

# Extract the data from the response
data = response.json()

# Extract the desired data
entries = data['data']['entries']

# Initialize an empty list to hold the entry data
entry_data = []

for entry in entries:
    entry_data.append({
        'ID': entry['rcsb_id'],
        'Initial Release Date': entry['rcsb_accession_info']['initial_release_date'],
        'PubMed ID': entry['rcsb_primary_citation']['pdbx_database_id_PubMed'],
        'DOI': entry['rcsb_primary_citation']['pdbx_database_id_DOI'],
        'Resolution': entry['rcsb_entry_info']['resolution_combined']
    })

# Create a DataFrame from the entry data
df = pd.DataFrame(entry_data)

# Write the DataFrame to a CSV file
df.to_csv('pdb_data_oprm1.csv', index=False)
