import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from collections import Counter


def load_data_from_json(file_path):
    """
    Loads data from a JSON file where objects are separated by newline.

    :param file_path: Path of the input file in JSON format
    :return: List of dictionaries
    """
    data = []
    with open(file_path, 'r') as file:
        print(file_path, 'successfully open')
        for line in file:
            line = line.strip()
            if line:
                data_line = json.loads(line)
                data.append(data_line)
    return data


def count_trait_gene_combinations(data):
    """
    Counts unique combinations of Trait ID and gene ID.

    :param data: List of dictionaries
    :return: Dictionary with Trait IDs as keys and counts of genes as values
    """
    traits_gene_count = {}
    unique_genes = set()
    for item in data:
        trait_id = item.get("trait ID", "")
        if trait_id not in traits_gene_count:
            unique_genes = set()
        gene_info = item.get("gene_info", [])
        for gene in gene_info:
            gene_id = [gene.get("gene_id")]
            unique_genes.update(gene_id)

        traits_gene_count[trait_id] = len(unique_genes)
    print(traits_gene_count)
    return traits_gene_count


def create_barplot(trait_gene_dict, th):
    """
    Creates a bar plot of counts of genes per Trait ID with LOD scores above a
    threshold of 3.

    :param trait_gene_dict: Dictionary with Trait IDs as keys and counts of
    genes as values
    """
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Gill Sans MT']
    filtered_data = [(trait, count) for trait, count in trait_gene_dict.items() if count > 3]

    # Check if filtered_data is empty
    if not filtered_data:
        return

    if th == 0:
        th_str = '95%'
    else:
        th_str = str(th)

    traits, counts = zip(*filtered_data)
    x_pos = range(len(traits))
    plt.bar(x_pos, counts, align='center', color='#9BAFB5',
            edgecolor='#536872', label='Genes Count')
    plt.xticks(x_pos, [trait.capitalize() for trait in traits], rotation=90)
    plt.xlabel('Trait ID')
    plt.ylabel(f'Gene Count (LOD > {th_str})')
    plt.title(f'Gene Count with LOD > {th_str} per Trait ID')
    plt.axhline(np.mean(counts), color='#475B61', linestyle='--',
                linewidth=2, alpha=0.6,
                label=f'Mean: {np.mean(counts):.2f}')
    plt.legend()
    plt.tight_layout()
    plt.show()
    #plt.savefig(r'qtl-analysis/qtl')


def count_unique_genes(data):
    """
    Counts the number of unique genes for each compound.

    :param data: Data loaded from JSON file
    :return: Dictionary with compound names as keys and counts of unique genes as values
    """
    unique_gene_counts = {}
    for compound, pathways in data.items():
        unique_genes = set()
        for pathway, genes in pathways.items():
            if not pathway.startswith('path:map011') and not pathway.startswith('path:map012'):
                unique_genes.update(genes)
        unique_gene_counts[compound] = len(unique_genes)
    return unique_gene_counts


def filter_compounds_with_genes(trait_gene_counts, unique_gene_counts_filtered):
    """
    Filters compounds where none of the counts is zero.

    :param trait_gene_counts: Dictionary with Trait IDs as keys and counts of genes as values
    :param unique_gene_counts_filtered: Filtered dictionary with compound names as keys and counts of unique genes as values
    :return: List of compounds
    """
    return {compound: count for compound, count in unique_gene_counts_filtered.items() if compound in trait_gene_counts}


def plot_unique_gene_counts(unique_gene_counts_filtered):
    """
    Plots the counts of unique genes for each compound.

    :param unique_gene_counts_filtered: Filtered dictionary with compound names as keys and counts of unique genes as values
    """
    if not unique_gene_counts_filtered:
        return

    mean_value = np.mean(list(unique_gene_counts_filtered.values()))
    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    x_pos = list(unique_gene_counts_filtered.keys())
    plt.bar(x_pos, unique_gene_counts_filtered.values(), color='#C96731', edgecolor='#974D25', linewidth=1.2)
    plt.xticks(x_pos, [trait.capitalize() for trait in x_pos], rotation=90)
    plt.axhline(mean_value, color='red', linestyle='--', alpha=0.3, label=f'Mean: {mean_value:.2f}')
    plt.legend()
    plt.title('Unique Gene Counts for Each Compound')
    plt.xlabel('Compound')
    plt.ylabel('Unique Gene Count')
    plt.tight_layout()
    plt.show()


def plot_gene_counts_per_compound(trait_gene_counts, unique_gene_counts_filtered, compounds, Title, fig_name):
    """
    Plots gene counts per compound.

    :param trait_gene_counts: Dictionary with Trait IDs as keys and counts of genes as values
    :param unique_gene_counts_filtered: Filtered dictionary with compound names as keys and counts of unique genes as values
    :param compounds: List of compounds to plot
    """
    bar_width = 0.35
    positions = np.arange(len(compounds))

    # Sort compounds based on the count of unique genes in metabolic pathways
    compounds_sorted = sorted(compounds, key=lambda x: unique_gene_counts_filtered[x])

    # Retrieve corresponding values for sorted compounds
    mQTL_gene_counts_values = [trait_gene_counts[compound] for compound in compounds_sorted]
    unique_gene_counts_filtered_values = [unique_gene_counts_filtered[compound] for compound in compounds_sorted]

    trait_gene_counts_mean = np.mean(mQTL_gene_counts_values)
    unique_gene_counts_filtered_mean = np.mean(unique_gene_counts_filtered_values)

    plt.figure(figsize=(12, 6))

    plt.bar(positions - bar_width / 2, mQTL_gene_counts_values, bar_width, color='#9BAFB5', edgecolor='#536872',
            linewidth=1.2, label='mQTL above thr.')

    plt.bar(positions + bar_width / 2, unique_gene_counts_filtered_values, bar_width, color='#C96731',
            edgecolor='#974D25', linewidth=1.2, label='Metabolic pathway')

    plt.axhline(trait_gene_counts_mean, color='blue', linestyle='--', alpha=0.5,
                label=f'Mean count (mQTL above thr.): {trait_gene_counts_mean:.2f}')

    plt.axhline(unique_gene_counts_filtered_mean, color='red', linestyle='--', alpha=0.5,
                label=f'Mean count (Metabolic pathway): {unique_gene_counts_filtered_mean:.2f}')

    plt.xlabel('Compounds',fontsize=16)
    plt.ylabel('Count of Genes',fontsize=16)
    plt.title(Title)
    plt.xticks(positions, compounds_sorted, rotation=90, ha='right',
               fontsize=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(fig_name + '.png', dpi=300)
    plt.show()

def main(study_name, th):
    name = study_name.split('_')[0]
    # Load data from JSON file
    #file_region = (r'Intermediate results/Serin_etal_2017_pheno_genes_in_region.json')
    directory = f'{study_name}_{th}'
    file_region = f'{study_name}_genes_in_region_th_{th}.json'

    outfile = f'{name}_{th}_general_count.json'

    figure_name = outfile[:-5] + '.png'
    if th == 0:
        th_str = '95%'
    else:
        th_str = str(th)

    title_fig = (f'Gene Count Comparison per Compound\n({name} dataset, '
                 f'LOD threshold = {th_str})')

    file_kegg = 'Genes_in_path_filtered_test_version_v2.json'
    data_region = load_data_from_json(fr"{directory}\{file_region}")

    # Open the file of gene/pathways and load its contents as a dictionary
    with open(file_kegg, 'r') as file:
        data_kegg = json.load(file)

    # Count unique combinations of Trait ID and genes
    trait_gene_qtl_counts = count_trait_gene_combinations(data_region)

    # Create and display the bar plot
    create_barplot(trait_gene_qtl_counts, th)

    # Count the number of unique genes for each compound
    unique_gene_paths_counts = count_unique_genes(data_kegg)

    # Filter out compounds with 0 unique genes
    unique_gene_paths_counts_filtered = filter_compounds_with_genes(trait_gene_qtl_counts, unique_gene_paths_counts)

    # Plot the counts of unique genes for each compound
    plot_unique_gene_counts(unique_gene_paths_counts_filtered)

    # Plot gene counts per compound
    compounds = list(unique_gene_paths_counts_filtered.keys())
    plot_gene_counts_per_compound(trait_gene_qtl_counts, unique_gene_paths_counts_filtered, compounds, title_fig, figure_name)

    # Save DataFrame to JSON file
    trait_gene_counts_values = [trait_gene_qtl_counts[compound] for compound in compounds]
    unique_gene_counts_filtered_values = [unique_gene_paths_counts_filtered[compound] for compound in compounds]
    df = pd.DataFrame(
        {'mQTL Gene Counts': trait_gene_counts_values, 'map Gene Counts': unique_gene_counts_filtered_values},
        index=compounds)
    df.to_json(rf"{directory}\{outfile}", orient='index')

if __name__ == '__main__':
    main(study_name='Joosen_etal_2013', th=0)
    main(study_name='Serin_etal_2017_pheno', th=0)