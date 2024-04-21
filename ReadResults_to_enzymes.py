"""
Created on 22/02/2024
This script reads the files from a specific folder.
The files are in CSV format and contain the column names:
    1. metabolic_path
    2. chromosome
    3. marker_ID
    4. gene_ID

The output of this script is a json file with a dictionary:

    {   Compound_name (str, from the file name) :

            {   metabolic_path (str) :
                list of enzymes (str) related to the genes in the pathway
            }
    }

"""

import os
import requests
import json
import time


def get_enzymes_for_gene(gene_id):
    """
    Get all the enzymes related to a gene
    :param gene_id: str
    :return: list of str
    """
    # Construct the URL for the KEGG REST API to retrieve pathway information for the gene
    url = f'http://rest.kegg.jp/link/enzyme/{gene_id}'

    # Make a GET request to the KEGG API
    response = requests.get(url)
    time.sleep(1)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response text to extract enzyme IDs
        enzyme_ids = [line.split('\t')[1].split(':')[-1] for line in response.text.split('\n') if line.strip()]
        return enzyme_ids
    else:
        print(f"Failed to retrieve enzyme information for gene: {gene_id}")
        return None


def main(study_name, th):

    th_str = str(th).replace('.', '_')
    # Create directory name
    folder_path = f'{study_name}_{th_str}'
    entries = os.listdir(folder_path)
    print(entries)


    # Check if directory already exists
    if not os.path.exists(folder_path):
        # If directory doesn't exist, create it
        os.mkdir(folder_path)

    outfile_name = f'result_enzymes_{study_name}_{th_str}.txt'
    outfile = open(f'{folder_path}/{outfile_name}', 'w')

    compounds = {}

    total = len(entries)
    idx = 0
    for file in entries:
        idx+=1
        print(f'{idx}/{total}')
        # Avoid doocuments in the folder that are other than csv
        if not file.endswith('.csv'):
            continue
        # fom the name of the file, extract the compound name
        compound = file.split('_')[-3]
        #print(compound)
        enzymes = {}

        infile = f'{folder_path}/{file}'
        #print(infile, os.path.exists(folder_path))
        with open(infile, 'r') as csv_file:

            for row in csv_file:
                #print(row)
                row = row.split()
                # Skip the first row (names of the columns)
                if row[3] == 'gene_ID':
                    continue

                path_id = row[0]
                gene_id = row[3]
                enzyme_ids = get_enzymes_for_gene(gene_id)
                # Avoid Global and Overview Pathways
                if path_id.startswith('path:map011') or path_id.startswith('path:map012'):
                    continue
                if path_id in enzymes.keys():
                    # Avoid duplications
                    for enzyme_id in enzyme_ids:
                        if enzyme_id not in enzymes[path_id]:
                            enzymes[path_id] += [enzyme_id]
                        else:
                            print(f"Duplicate enzyme for {path_id}, {enzyme_id}")

                else:
                    enzymes[path_id] = enzyme_ids

        compounds[compound] = enzymes

    json.dump(compounds, outfile, indent=4)

if __name__ == '__main__':
    th = 0
    study_name = 'Joosen_etal_2013'
    main(study_name, th)
    study_name = 'Serin_etal_2017_pheno'
    main(study_name, th)




