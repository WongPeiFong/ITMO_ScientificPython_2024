import requests
import subprocess
import json
from Bio import SeqIO
import re

UNIPROT_BASE_URL = 'https://rest.uniprot.org/uniprotkb/accessions'
ENSEMBL_BASE_URL = 'https://rest.ensembl.org/sequence/id/'

UNIPROT_PATTERN = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]\d[A-Z][A-Z\d]{2}\d{1,2}')
ENSEMBL_PATTERN = re.compile(r'ENS[A-Z]{1,2}[A-Z]{3}\d{11}')

def get_uniprot(ids):
    accessions = ','.join(ids)
    http_args = {'params': {'accessions': accessions}}
    response = requests.get(UNIPROT_BASE_URL, **http_args)
    return response

def parse_response_uniprot(response):
    data = response.json()
    parsed_data = []
    for entry in data.get('results', []):
        parsed_data.append({
            'Entry': entry['primaryAccession'],
            'Protein names': entry['proteinDescription']['recommendedName']['fullName']['value'] if 'recommendedName' in entry['proteinDescription'] else None,
            'Gene names': ', '.join([gene['geneName']['value'] for gene in entry['genes']]) if 'genes' in entry else None,
            'Organism': entry['organism']['scientificName'],
            'Length': entry['sequence']['length']
        })
    return parsed_data

def get_ensembl(ids):
    headers = {'Content-Type': 'application/json'}
    data = {'ids': ids}
    response = requests.post(ENSEMBL_BASE_URL, headers=headers, data=json.dumps(data))
    return response

def parse_response_ensembl(response):
    data = response.json()
    parsed_data = {entry['id']: entry for entry in data}
    return parsed_data

def fetch_and_parse_data(ids_dict):
    uniprot_ids = ids_dict.get('Uniprot', [])
    ensembl_ids = ids_dict.get('ENSEMBL', [])

    parsed_data = {}
    if uniprot_ids:
        uniprot_response = get_uniprot(uniprot_ids)
        if uniprot_response.status_code == 200:
            parsed_data['Uniprot'] = parse_response_uniprot(uniprot_response)
        else:
            parsed_data['Uniprot'] = f"Failed to fetch data from Uniprot. Status code: {uniprot_response.status_code}"

    if ensembl_ids:
        ensembl_response = get_ensembl(ensembl_ids)
        if ensembl_response.status_code == 200:
            parsed_data['ENSEMBL'] = parse_response_ensembl(ensembl_response)
        else:
            parsed_data['ENSEMBL'] = f"Failed to fetch data from ENSEMBL. Status code: {ensembl_response.status_code}"

    return parsed_data

def find_uniprot_ensembl_ids(sequences, file_type):
    ids = {'Uniprot': [], 'ENSEMBL': []}
    pattern = UNIPROT_PATTERN if file_type == 'Protein' else ENSEMBL_PATTERN

    for sequence in sequences:
        match = re.search(pattern, sequence['description'])
        if match:
            ids_key = 'Uniprot' if file_type == 'Protein' else 'ENSEMBL'
            ids[ids_key].append(match.group())

    return ids

def call_seqkit_stats(fasta_file):
    try:
        result = subprocess.run(["seqkit", "stats", fasta_file], capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return str(e)

def parse_seqkit_stats(stats_output):
    if "DNA" in stats_output:
        return "DNA"
    elif "Protein" in stats_output:
        return "Protein"
    else:
        return None 

def process_fasta_file(fasta_file):
    stats_output = call_seqkit_stats(fasta_file)
    if "Error" in stats_output:
        return stats_output

    file_type = parse_seqkit_stats(stats_output)
    if file_type is None:
        return "Unable to determine file type."

    sequences = []
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append({
                "description": record.description,
                "sequence": str(record.seq)
            })

    ids = find_uniprot_ensembl_ids(sequences, file_type)
    
    if ids['Uniprot'] or ids['ENSEMBL']:
        db_info = fetch_and_parse_data(ids)
    else:
        db_info = "No Uniprot or ENSEMBL IDs found in sequence descriptions."

    final_output = {
        "seqkit_stats": stats_output,
        "sequences_info": sequences,
        "database_info": db_info,
        "database_name": "Uniprot" if file_type == 'Protein' else 'ENSEMBL'
    }
    return final_output

fasta_file = "hw_file1.fasta"
result = process_fasta_file(fasta_file)
print(result)
