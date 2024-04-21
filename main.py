"""
This script processes QTL data from a given study and performs gene pathway
analysis.

Parameters:
    study_name (str): The name of the study containing QTL data.

Returns:
    None

Dependencies:
    - matplotlib.pyplot as plt
    - pandas as pd
    - json
    - compare_genes (custom module)
    - QTL (custom module)
    - os
    - time

Outputs:
    - CSV files containing comparison results between genes and pathways.
    - Plots of LOD scores with markers above a certain threshold.

Example:
    To run this script, simply call the `main` function with the name of the study:
        main("study_name")

    Ensure that all necessary files and modules are in the working directory.
"""

import matplotlib.pyplot as plt
import pandas as pd
import json
import compare_genes
import QTL
import os
import time

def main(study_name, th):
    # Set initial patameters:
    #th = 3.85
    th_str = str(th).replace('.', '_')
    # Create directory name
    folder = f'{study_name}_{th_str}'

    # Check if directory already exists
    if not os.path.exists(folder):
        # If directory doesn't exist, create it
        os.mkdir(folder)

    # Get the QTLs from a Study (Object)
    study = QTL.Study(study_name)

    # get the ID's of the traits as used in the study
    traits_codes = study.LOD_scores.index
    # Open the list of traits of interest, so we avoid the traits that are not
    trait_df = pd.read_csv('List of traits.csv',delimiter=';',header=None)

    # Open the document with the dictionary of trait : {paths: genes}
    # This document contains all the pathways in KEGG related to a trait,
    # And all the genes related to that pathway
    with open('Genes_in_path_filtered_test_version_v2.json', 'r') as genes_in_path:

        # reconstruct the str data to a dictionary type
        Genes_in_path = json.load(genes_in_path)

    # Iterate over the traits from the study
    total = len(list(traits_codes))
    i = 1   # Save an index in case we need to pause the run

    for trait in traits_codes:
        print (f'{i}/{total}')
        # We will compare the name of the trait to make sure the trait is
        # in the list of traits of interest
        trait_row = trait_df[trait_df[0] == trait]
        # conditional statement
        # to see if this trait is associated with a trait ID or not
        if not trait_row.empty:
            #print(trait_row)
            trait_code = trait_row[0].values[0]
            #print('trait = ', trait_code)
            trait_ID = trait_row[1].values[0]
            # conditional to see if the trait_ID is present in Genes_in_path
            # some traits are not in KEGG, so we can not use them
            # so, we will contunue only with the ones that are useful
            if trait_ID in Genes_in_path:
                # Get pathways related to trait (list)
                genes_path = Genes_in_path[trait_ID]
                # get the LOD and marker of a peak associated with the trait_code
                peak_LOD, peak_markers = study.get_peak(
                    trait_code,
                    threshold=th
                )

                file = f'{study_name}_{i}_{trait_ID}_{th_str}_comparison'
                # The function compare_genes writes the results to a CSV file.
                file_path = f'{folder}\\{file}'
                resulting_genes, high_markers = compare_genes.compare_genes(
                    genes_path,
                    peak_markers,
                    file_path
                )

                study.plot_scores(
                    trait_code,
                    marker=high_markers,
                    threshold=th,
                    fig_title=f'{trait_ID}'
                )

                figure = f'{study_name[:3]}_{i}_{trait_ID}_{th_str}'
                path = f'{study_name}\\{figure}.png'
                plt.savefig(path, dpi=300, bbox_inches='tight')
                time.sleep(1)


                plt.close()
        i += 1



if __name__ == '__main__':
    main('Serin_etal_2017_pheno', 0)
    #main('Joosen_etal_2013', 0)

