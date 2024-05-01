import requests
import subprocess
from Bio import SeqIO
import re

def get_uniprot(ids):
    base_url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': 'ID',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(ids)
    }
    response = requests.get(base_url, params=params)
    return response

def parse_response_uniprot(response):
    parsed_data = []
    lines = response.text.split('\n')
    for line in lines[1:]:
        fields = line.strip().split('\t')
        parsed_data.append({
            'Entry': fields[0],
            'Status': fields[1],
            'Entry name': fields[2],
            'Protein names': fields[3],
            'Gene names': fields[4],
            'Organism': fields[5],
            'Length': fields[6]
        })
    return parsed_data

def get_ensembl(ids):
    base_url = 'https://rest.ensembl.org/sequence/id/'
    headers = {'Content-Type': 'application/json'}
    data = {'ids': ids, 'content-type': 'application/json'}
    response = requests.post(base_url, headers=headers, data=json.dumps(data))
    return response

def parse_response_ensembl(response):
    parsed_data = {}
    data = response.json()
    for entry in data:
        parsed_data[entry['id']] = {
            'sequence': entry['seq'],
            'molecule': entry['molecule'],
            'description': entry['desc']
        }
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
uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]\d[A-Z][A-Z\d]{2}\d{1,2}$')
ensembl_pattern = re.compile(r'^ENS[A-Z]{1,2}[A-Z]{3}\d{11}$')

def find_uniprot_ensembl_ids(sequences, file_type):
    ids = {'Uniprot': [], 'ENSEMBL': []}
    pattern = uniprot_pattern if file_type == 'Protein' else ensembl_pattern

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

def find_uniprot_ensembl_ids(description, file_type):
    uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]\d[A-Z][A-Z\d]{2}\d{1,2}$')
    ensembl_pattern = re.compile(r'^ENS[A-Z]{1,2}[A-Z]{3}\d{11}$')

    pattern = uniprot_pattern if file_type == 'Protein' else ensembl_pattern
    match = re.search(pattern, description)

    if match:
        return match.group()
    else:
        return None

def call_database_api(ids, database):
    if database == 'Uniprot':
        base_url = 'https://www.uniprot.org/uploadlists/'
        params = {
            'from': 'ID',
            'to': 'ACC',
            'format': 'tab',
            'query': ' '.join(ids)
        }
        response = requests.get(base_url, params=params)
        return response.text.strip()
    elif database == 'ENSEMBL':
        base_url = 'https://rest.ensembl.org/sequence/id/'
        headers = {'Content-Type': 'application/json'}
        data = {'ids': ids, 'content-type': 'application/json'}
        response = requests.post(base_url, headers=headers, data=json.dumps(data))
        return response.json()
    else:
        return f"Invalid database: {database}"

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

    ids = []
    for sequence in sequences:
        id_match = find_uniprot_ensembl_ids(sequence['description'], file_type)
        if id_match:
            ids.append(id_match)

    if ids:
        db_info = call_database_api(ids, 'Uniprot' if file_type == 'Protein' else 'ENSEMBL')
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
