#  first task
import requests
import json
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
uniprot_ids = ['P11473', 'Q91XI3']
ensembl_ids = ['ENSMUSG00000041147', 'ENSG00000139618']
uniprot_response = get_uniprot(uniprot_ids)
parsed_uniprot_data = parse_response_uniprot(uniprot_response)
print("Uniprot Data:")
print(parsed_uniprot_data)

ensembl_response = get_ensembl(ensembl_ids)
parsed_ensembl_data = parse_response_ensembl(ensembl_response)
print("\nENSEMBL Data:")
print(parsed_ensembl_data)

#  second task
import re
uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9][A-Z][A-Z\d]\d$')
ensembl_pattern = re.compile(r'^ENS[A-Z]{3}\d{11}$')
def fetch_and_parse_data(ids):
    uniprot_ids = []
    ensembl_ids = []
    for id in ids:
        if uniprot_pattern.match(id):
            uniprot_ids.append(id)
        elif ensembl_pattern.match(id):
            ensembl_ids.append(id)
        else:
            print(f"ID '{id}' does not match any database pattern.")

    parsed_data = {}
    if uniprot_ids:
        uniprot_response = get_uniprot(uniprot_ids)
        parsed_uniprot_data = parse_response_uniprot(uniprot_response)
        parsed_data['Uniprot'] = parsed_uniprot_data
    if ensembl_ids:
        ensembl_response = get_ensembl(ensembl_ids)
        parsed_ensembl_data = parse_response_ensembl(ensembl_response)
        parsed_data['ENSEMBL'] = parsed_ensembl_data
    return parsed_data
ids = ['ENSMUSG00000041147', 'ENSG00000139618', 'P11473', 'Q91XI3']
parsed_data = fetch_and_parse_data(ids)
print("Parsed Data:")
print(parsed_data)
