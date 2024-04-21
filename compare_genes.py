import GenesInRegion
import QTL
import pandas as pd


def compare_genes(dict_path_genes, df_markers, out_file):
    """
    Compare genes in a region with predefined pathways.
    The function writes the results to a CSV file.

    Parameters:
    - dict_path_genes: Dictionary of pathways and associated gene lists.
                        {'path:mapXXXXX': [list of genes]}
    - df_markers: Pandas DataFrame containing marker information.
    - out_file = str, name for a new file to retrieve the returns

    Returns:
    - result_list: List of dictionaries with information about found genes,
            associated pathways, chromosome, region start and end, and marker.
    - file .csv that contains metabolic path, chromosome, marker_ID and gene_ID
    """
    # initialize
    result_dictionary = {}
    markers = []
    out_file = out_file + '.csv'
    file = open(out_file, 'w')
    file.write('metabolic_path\tchromosome\tmarker_ID\tgene_ID\n')

    GIR = GenesInRegion.GenesInRegion()
    if type(df_markers) == pd.DataFrame:
        totalN = len(df_markers.index)
        if totalN < 1:
            return result_dictionary, set(markers)
        idx = 1
        for n in df_markers.index:
            marker_ID = n #Serin dataset
            #marker_ID = df_markers.loc[n].marker #Joosen dataset
            print(marker_ID)

            chrom, start_p, end_p = QTL.get_region(df_markers, n)
            genes_reg = GIR.genes_in_region(chrom, start_p, end_p)
            #GIR.daysTillDeadline = 9
            print('marker {0}/{1}'.format(idx,totalN))


            for gene in genes_reg:
                gene_id = gene['gene_id']
                # Check if the value is present in any of the lists within the
                # dictionary values
                for path in dict_path_genes:
                    gene_list = dict_path_genes[path]
                    if gene_id in gene_list:
                        file.write(f'{path}\t{chrom}\t{marker_ID}\t{gene_id}\n')
                        result_dictionary[gene_id] = (path, chrom, marker_ID)
                        markers += [str(marker_ID)]
            idx += 1
        file.close()

    return result_dictionary, set(markers)
