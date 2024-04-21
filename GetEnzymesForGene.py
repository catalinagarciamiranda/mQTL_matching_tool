"""
This Python module contains a function `get_enzymes_for_gene` that retrieves
a list of enzyme IDs related to a given gene ID from the KEGG REST API.

This module requires the `requests` library to make HTTP requests.

Functions:
    get_enzymes_for_gene(gene_id: str) -> List[str] or None
"""

import requests

def get_enzymes_for_gene(gene_id):
    """
    This function retrieves all the enzyme IDs related to a given gene ID by
    making a GET request to the KEGG REST API.

    :param gene_id: str, for example 'ath:AT1G06290', gene ID of Ath.
    :return: list[str], for example ['ec:1.3.3.6'], all the enzyme IDs
    related to the gene ID. If the request fails, it returns None.

    The function constructs the URL for the KEGG REST API to retrieve pathway
    information for the gene. It then makes a GET request to the KEGG API.
    If the request is successful (status code 200), it parses the response text
    to extract enzyme IDs. The enzyme IDs are extracted from the response text
    by splitting each line on the tab character ('\t'), taking the second
    element, splitting this element on the colon character (':'), and taking
    the last element. If the request is not successful, it prints an error
    message and returns None.
    """
    # Construct the URL for the KEGG REST API to retrieve pathway information
    # for the gene
    url = f'http://rest.kegg.jp/link/enzyme/{gene_id}'

    # Make a GET request to the KEGG API
    response = requests.get(url)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response text to extract enzyme IDs
        enzyme_ids = [line.split('\t')[1].split(':')[-1]
                      for line in response.text.split('\n') if line.strip()]
        return enzyme_ids

    else:
        print(f"Failed to retrieve enzyme information for gene: {gene_id}")
        return None


# Test
if __name__ == "__main__":
    # initializes an empty dictionary
    enzymes = {}
    # opens the CSV file in read mode
    file_name = 'Serin_etal_2017_pheno_3_9/' \
                'Serin_etal_2017_pheno_10_Citric acid_3_9_comparison.csv'
    csv_file = open(file_name, 'r')

    # start a loop that iterates over each row in the CSV file
    for row in csv_file:
        # split the current row into a list of strings
        row = row.split()

        # avoid the first line
        if row[3] == 'gene_ID':
            continue

        # assign the corresponding elements to the variables path and gene_id
        path = row[0]
        gene_id = row[3]
        # call the function get_enzymes_for_gene
        enzyme_ids = get_enzymes_for_gene(gene_id)

        # check if path is already a key in the enzymes dictionary
        # If it is, it appends enzyme_ids to the existing list of enzyme IDs
        # for that path
        if path in enzymes.keys():
            enzymes[path] += enzyme_ids
        # If path is not already a key in the enzymes dictionary, this line
        # creates a new entry in the dictionary with path as the key and
        # enzyme_ids as the value.
        else:
            enzymes[path] = enzyme_ids

    # print an index and the enzyme IDs for each path
    n = 1
    for element in enzymes.values():
        print(n, element)
        n += 1
