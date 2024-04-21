"""
Updated on 04/04/2024
@author: Catalina García Miranda

This script automates the retrieval of pathway IDs from the KEGG pathway
database for a list of metabolites.
It handles potential errors such as failed API requests and retries failed
requests until successful.
The retrieved pathway data is stored in a JSON file for further analysis.

Functions:
    - get_kegg_pathways(compound_name): Retrieves pathways related to a given
    compound from KEGG.
    - append_to_json(new_data, file_path): Appends new data to a JSON file.

Main Execution:
    - Reads a list of metabolites from a CSV file.
    - Retrieves pathways for each metabolite.
    - Saves the data to a JSON file.
    - Retries failed requests for compounds until successful.

"""

import pandas as pd
import requests
import json
import logging
import urllib.parse
import time
import re

def get_kegg_pathways(compound_name):
    """
    Retrieves a list of pathway IDs from the KEGG pathway database related to a
    given compound name.

    :param compound_name: str
        A metabolite, for instance 'Raffinose'.
    :return: list of str
        List of pathway IDs related to the compound_name according to KEGG API.
        For example ['path:map00052', 'path:map02010'].
    """

    # Set the URL of KEGG
    base_url = 'http://rest.kegg.jp/'
    # Replace spaces with '%20' in compound_name
    encoded_compound_name = urllib.parse.quote(compound_name)
    search_url = f"{base_url}find/compound/{encoded_compound_name}"
    #print(search_url)

    # Search for the compound
    response = requests.get(search_url)

    # Handle the case where the server did not respond successfully

    if response.status_code != 200:
        t = 0
        while t < 3 and response.status_code != 200:
            t += 1
            response = requests.get(search_url)

    if response.status_code != 200:
        return []

    # Extract lines of text from the response
    lines = response.text.strip().split("\n")

    if not lines:
        #logging.info("No compounds found for %s.", compound_name)
        return []

    pattern = pattern = re.compile(rf'(?<![\w-])(?:[LD]-)?{compound_name}\b(?!-)', flags=re.IGNORECASE)
    #print(pattern)
    # Initialize an empty list to save the pathways
    pathways = []
    # search the compound name in the response text
    for line in lines:
        compound_id, compound_names = line.split("\t")
        # I want to keep a compound ID only if it is associated to exactly the
        # comound name, for instance, if the compound is "Alanine", I do not
        # want to keep "L-gamma-Glutamyl-D-alanine"
        # \bAlanine\b matches the word "Alanine" as a whole word,
        # ensuring that it is not part of a longer word (e.g., "D-Alanine")
        finder = re.search(pattern, compound_names)
        if finder:
            #print(compound_id)
            # Search list of pathways related to the compound ID
            pathway_url = f"{base_url}link/pathway/{compound_id}"
            #print(pathway_url)

            response = requests.get(pathway_url)
            # Handle the case where the server did not respond successfully
            if response.status_code != 200:
                t = 0
                while t < 3 and response.status_code != 200:
                    t += 1
                    response = requests.get(pathway_url)

            if response.status_code != 200:
                return set(pathways)

            # Extract lines of text from the response
            resp_lines = response.text.strip().split("\n")

            if lines:
                for r_line in resp_lines:
                    if "\t" in r_line:
                        pathway_id = r_line.split("\t")[1]

                        # Filter out Global and Overview pathways
                        if not pathway_id.startswith('path:map011') and not pathway_id.startswith('path:map012'):
                            pathways.append(pathway_id)
                        #else:
                            #logging.debug("Filtered pathway ID: %s", pathway_id)

    return set(pathways)  # Convert to set to remove duplicates

def append_to_json(new_data, file_path):
    """Append new data to existing JSON file."""
    with open(file_path, 'r+') as file:
        try:
            data = json.load(file)
        except json.JSONDecodeError:
            data = {}  # File is empty or not valid JSON
        data.update(new_data)
        file.seek(0)
        json.dump(data, file, indent=4)


if __name__ == "__main__":

    # Configure logging
    logging.basicConfig(level=logging.INFO)

    # Read a list of metabolites from a CSV file
    file_name = 'List of traits.csv'
    trait_df = pd.read_csv(file_name, delimiter=';', header=None)

    # Create a text document to save the data
    json_file = 'pathways_filtered_test_version_v2.json'

    # Initialize a dictionary to save the data
    path_dict = {}

    # Save the dictionary to a text file in JSON format
    json.dump(path_dict, open(json_file, 'w'), indent=4)

    failed_compounds = set()
    consecutive_failures = 0

    # Initial flag to run the function for all traits
    run_for_all_traits = True
    # Loop until either 3 consecutive failures occur or the failure set remains unchanged for two iterations
    while consecutive_failures <= 3:
        # Run for all traits or only for failed compounds
        if run_for_all_traits:
            compounds_to_process = set(trait_df.iloc[:, 1])
        else:
            compounds_to_process = set(failed_compounds)

        failed_compounds = set()
        # Iterate over the compounds
        for kegg_comp in compounds_to_process:
            # Find all the pathways related to a metabolite
            #logging.info("Searching for pathways related to metabolite: %s", kegg_comp)
            pathways = get_kegg_pathways(str(kegg_comp))

            # Check if pathways retrieval was successful
            if pathways:
                path_dict[str(kegg_comp)] = list(pathways)
            else:
                #logging.error("Failed to retrieve pathways for metabolite: %s", kegg_comp)
                failed_compounds.add(str(kegg_comp))  # Store the failed compound
                print("List of failed compounds:", failed_compounds)

        # increment consecutive_failures counter
        if failed_compounds == compounds_to_process:
            consecutive_failures += 1
        else:
            consecutive_failures = 0

        print(f'failures={consecutive_failures}/3')

        # Change the flag to run for failed compounds only in the next iteration
        run_for_all_traits = False

        append_to_json(new_data=path_dict, file_path=json_file)

    # Print the list of failed compounds
    print("List of failed compounds:", failed_compounds)
