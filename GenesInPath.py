"""Updated version og GenesInPath 05/04/2024
@author: Catalina García Miranda

This script retrieves Arabidopsis gene IDs related to a list of KEGG pathways.

It reads a JSON file containing a dictionary where,
 - keys are metabolites, and
 - values are lists of KEGG pathway IDs relevant to each metabolite.
For each metabolite, it retrieves the gene IDs associated with the pathways
using the KEGG REST API.

The retrieved data is stored in another JSON file with the following structure:
{
    "metabolite1": {
        "pathway_ID1": ["gene_ID1", "gene_ID2", ...],
        "pathway_ID2": ["gene_ID3", "gene_ID4", ...],
        ...
    },
    "metabolite2": {
        ...
    },
    ...
}

Note: This script implements a retry mechanism for HTTP requests that do not
return a status code of 200.
"""

import requests
import time
import json


def get_response(url):
    """
    Retrieve a response from the specified URL, retrying a maximum of three
    times if the initial request fails.

    :param url: str
        The URL to send the request to.

    :return: requests.Response or None
        The response object if the request is successful (status code 200),
        otherwise None.
    """
    max_retries = 3
    retries = 0
    while retries < max_retries:
        response = requests.get(url)
        if response.status_code == 200:
            return response
        else:
            retries += 1
            print(f"Request failed with status code {response.status_code}. "
                  f"Retrying...")
            time.sleep(1)  # Add a small delay before retrying
    # If all retries fail, return None
    return None


def get_genes_for_pathways(pathway_list):
    """Retrieve Arabidopsis gene IDs related to a list of KEGG pathway IDs.

       :param: pathway_list (list of str): List of KEGG pathway IDs,
           e.g., ["path:map01061", "path:map01062"].

       :return: dict: A dictionary where,
            - keys are pathway IDs, and
            - values are lists of gene IDs associated with each pathway.
            - If no genes are found for a pathway, it is not included.
       """

    genes = {}
    for pathway_id in pathway_list:
        if (pathway_id.startswith('path:map011')
                or pathway_id.startswith('path:map012')):
            print('filtered pathway ID: ', pathway_id)
            continue
        pathway_ath = 'ath' + pathway_id[8:]
        pathway_url = f"https://rest.kegg.jp/link/ath/{pathway_ath}"

        # Use get_response function to handle retries
        response = get_response(pathway_url)
        if response is not None:
            text_data = response.text
            lines = text_data.strip().split('\n')
            genes[pathway_id] = []
            for line in lines:
                columns = line.split('\t')
                if len(columns) > 1:
                    second_column = columns[1].strip()
                    genes[pathway_id].append(second_column)
            if not genes[pathway_id]:
                del genes[pathway_id]
    return genes


def main():
    """Main function to retrieve gene IDs for pathways related to metabolites.

    1. Opens a JSON file containing a dictionary where,
    - keys are metabolites, and
    - values are lists of KEGG pathway IDs.
    2. Iterates through the metabolites, retrieves gene IDs associated with
    their pathways, and stores the data in another JSON file.

    Note:
        This function implements a retry mechanism for HTTP requests that do
        not return a status code of 200.
    """
    path_file = 'pathways_filtered_test_version_v2.json'
    mtbl_dict = {}
    outfile = open('Genes_in_path_filtered_test_version_v2.json', 'w')
    with open(path_file, 'r') as path_file:
        dict_comp_path = json.load(path_file)
        for metbl in dict_comp_path:
            print(metbl)
            print('working on it')
            pathway_list = dict_comp_path[metbl]
            print('len pathway list = ', len(pathway_list))
            genes_dict = get_genes_for_pathways(pathway_list)
            mtbl_dict[metbl] = genes_dict
        json.dump(mtbl_dict, outfile, indent=4)
    outfile.close()


if __name__ == '__main__':
    main()
