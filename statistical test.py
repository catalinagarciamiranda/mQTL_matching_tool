import json
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def main(study_name, th):
    name = study_name.split('_')[0]
    print(name)
    if th == 0:
        th_str = '95%'
    else:
        th_str = str(th)

    title = (f'Number of Matching Genes per Compound in {name} '
             f'dataset (LOD threshold = {th_str}),\n'
             f'Significance Level (Fisher Exact p-value adj. < 0.05)')

    directory = f'{study_name}_{th}'
    # file to save the figure in
    file_figure = f'{name}_{th}_significant_genes.png'

    # Load gene counts data from gene_counts.json
    file_general_count = rf'{name}_{th}_general_count.json'

    with open(rf'{directory}\{file_general_count}', 'r') as file:
        general_counts_data = json.load(file)

    # Convert gene counts data to DataFrame
    gene_counts_df = pd.DataFrame(general_counts_data).T.reset_index()
    gene_counts_df.columns = ['Compound', 'mQTL Gene Counts', 'map Gene Counts']
    print(len(gene_counts_df['Compound']))

    # Load filtered gene counts data from filtered_resulting_gene_count.json
    file_filtered_matched = rf'{study_name}_{th}/{name}_{th}_count_matched.json'
    #print(file_filtered_matched)

    with open(file_filtered_matched, 'r') as file:
        filtered_gene_counts_data = json.load(file)

    # Convert filtered gene counts data to DataFrame
    filtered_gene_counts_df = pd.DataFrame(filtered_gene_counts_data)
    print(filtered_gene_counts_df)

    # Merge gene counts DataFrames
    # Merge gene counts DataFrames while keeping all rows and fill missing values with 0
    merged_df = pd.merge(gene_counts_df, filtered_gene_counts_df, on='Compound')
    merged_df.fillna(0, inplace=True)

    # Add 'Total Genes' column with fixed value
    merged_df['Total Genes'] = 27533 -231

    #print("Merged DataFrame:")
    print(merged_df)

    # Create an empty list to store the p-values
    p_values = []

    # Iterate through each row in the DataFrame
    for index, row in merged_df.iterrows():
        # Extract the relevant counts for the contingency table
        a = row['map Gene Counts']  # Number of genes in category c1 ('map Gene Counts')
        b = row['Total Genes'] - row['map Gene Counts']  # Number of genes not in category c1
        c = row['Gene_Count']  # Number of genes in both categories c1 and c2 ('Gene_Count')
        d = row['mQTL Gene Counts'] - row['Gene_Count']  # Number of genes in category c2 ('mQTL Gene Counts') but not in c1

        # Perform Fisher exact test
        oddsratio, p_value = fisher_exact([[a, b], [c, d]], alternative='two-sided')

        # Append the p-value to the list
        p_values.append(p_value)

    # Add the list of p-values as a new column in the DataFrame
    merged_df['Fisher_Exact_p_value'] = p_values


    # Compute adjusted p-values using Benjamini-Hochberg procedure
    reject, adjusted_p_values, _, _ = multipletests(merged_df['Fisher_Exact_p_value'], method='fdr_bh')

    # Add adjusted p-values to the DataFrame
    merged_df['Adj_p_value'] = adjusted_p_values

    # Print the DataFrame to verify the changes
    print(merged_df)

    ##############################################################################

    import matplotlib.pyplot as plt
    import seaborn as sns

    # Define significance threshold
    significance_threshold = 0.05  # You can adjust this value as needed

    # Filter the DataFrame for significant p-values
    significant_df = merged_df[merged_df['Adj_p_value'] < significance_threshold]

    # Sorted
    merged_df = merged_df.sort_values(by='Gene_Count', ascending=True)
    # Get all unique compounds
    all_compounds = merged_df['Compound'].unique()


    # Plotting
    plt.figure(figsize=(12, 6))

    # Plot all genes
    sns.barplot(data=merged_df, x='Compound', y='Gene_Count', color='#9BAFB5',
                label='Not significant', order=all_compounds)

    # Highlight significant genes
    sns.barplot(data=significant_df, x='Compound', y='Gene_Count', color='#F6A21D',
                label='Significant', order=all_compounds)
    #plt.title(title, fontsize=16)
    plt.xlabel('Compound', fontsize=16)
    plt.ylabel('Gene Count', fontsize=16)
    plt.grid(axis='y', alpha=0.5)
    # Capitalize the first character of each label on the x-axis
    plt.xticks(rotation=90, ticks=range(len(all_compounds)),
               labels=[compound.capitalize() for compound in all_compounds],
               fontsize=14)

    plt.legend(loc='upper left', fontsize=14)
    plt.ylim(0, 70)  # Set y-axis limits
    plt.tight_layout()  # Adjust layout to prevent overlap
    plt.savefig(file_figure)
    plt.show()

if __name__ == '__main__':
    main(study_name='Joosen_etal_2013', th=0)
    main(study_name='Serin_etal_2017_pheno', th=0)
