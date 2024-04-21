import gzip
import re
import QTL
import pandas as pd
import json

class GenesInRegion:
    def __init__(self):
        self.GFF_file = []
    def genes_in_region(self, chromosome, start_position, end_position):
        chromosome = f'Chr{chromosome}'
        # Define the path to the GFF file
        gff_file_path = r'C:\Users\Gebruiker\Documents\Estudio\Master Bioinformatics\Thesis\qtl-analysis\Araport11_GFF3_genes_transposons.Jul2023.gff.gz'

        # Define a regular expression pattern to match attribute key-value pairs
        attribute_pattern = re.compile(r'([^=;]+)=("[^"]+"|[^;]+)')

        # Initialize a list to store the genes in the specified range
        genes_in_range = []

        # Open and read the GFF file
        if not self.GFF_file:
            with gzip.open(gff_file_path, 'rt') as file:
                for line in file:
                    # Skip comments
                    if line.startswith('#'):
                        continue

                    # Split the line into columns
                    columns = line.strip().split('\t')
                    self.GFF_file.append(columns)

                    # Check if the line represents a gene on the chromosome of interest
                    # and within the specified range
                    if len(columns) > 2 \
                            and columns[2] == 'gene' \
                            and columns[0] == chromosome:
                        start = int(columns[3])
                        end = int(columns[4])
                        if start_position <= start <= end_position \
                                or start_position <= end <= end_position:
                            # Extract and parse the attributes to get
                            # gene_id and gene_name
                            attributes = columns[8]
                            attr_dict = {match.group(1): match.group(2).strip('"')
                                         for match
                                         in attribute_pattern.finditer(attributes)}
                            gene_id = 'ath:' + attr_dict.get('ID')
                            gene_name = attr_dict.get('Name')

                            # Add the gene information to the list
                            genes_in_range.append({
                                'gene_id': gene_id,
                                'gene_name': gene_name,
                                'chr': columns[0],
                                'start': start,
                                'end': end
                            })
        else:
            # Check if the line represents a gene on the chromosome of interest
            # and within the specified range
            for columns in self.GFF_file:
                if len(columns) > 2 and columns[2] == 'gene' and columns[0] == chromosome:
                    start = int(columns[3])
                    end = int(columns[4])
                    if start_position <= start <= end_position or start_position <= end <= end_position:
                        # Extract and parse the attributes to get
                        # gene_id and gene_name
                        attributes = columns[8]
                        attr_dict = {match.group(1): match.group(2).strip('"') for match in
                                     attribute_pattern.finditer(attributes)}
                        gene_id = 'ath:' + attr_dict.get('ID')
                        gene_name = attr_dict.get('Name')

                        # Add the gene information to the list
                        genes_in_range.append(
                            {'gene_id': gene_id, 'gene_name': gene_name, 'chr': columns[0], 'start': start, 'end': end})

        # Print the genes found in the specified range

        # print(f'Number of genes in {chromosome} between {start_position} and {end_position}: ', len(genes_in_range))

        return genes_in_range

def main(study_name, th):
    th_str = str(th).replace('.', '_')
    file = f'{study_name}_genes_in_region_th_{th_str}.json'
    study = QTL.Study(study_name)
    traits_codes = study.LOD_scores.index
    total = len(traits_codes)
    trait_df = pd.read_csv(r'C:\Users\Gebruiker\Documents\Estudio\Master Bioinformatics\Thesis\qtl-analysis\List of traits.csv', delimiter=';', header=None)

    with open(file, 'w') as outfile:
        # json.dump(('trait ID', ['gene_id', 'gene_name', 'chr', 'start', 'end']), outfile)

        i = 0
        for trait in traits_codes:
            i += 1
            print(f'{i}/{total}')
            trait_row = trait_df[trait_df[0] == trait]

            if not trait_row.empty:
                trait_code = trait_row[0].values[0]
                trait_ID = trait_row[1].values[0]

                peak_LOD, peak_markers = study.get_peak_local_max(trait_code, threshold=th)
                subtotal = len(peak_markers.index)
                j = 0
                for n in peak_markers.index:
                    j += 1
                    print(f'{j}/{subtotal}')
                    chr, start_p, end_p = QTL.get_region(peak_markers, n)
                    GIR = GenesInRegion()
                    genes_in_range = GIR.genes_in_region(chr, start_p, end_p)
                    #json.dump((trait_ID, genes_in_range), outfile)
                    json.dump({"trait ID": trait_ID, "gene_info": genes_in_range}, outfile)
                    outfile.write('\n')

if __name__ == '__main__':
    main(study_name='Joosen_etal_2013', th=0)
    #main(study_name='Joosen_etal_2013', th=3)
    #main(study_name='Serin_etal_2017_pheno', th=3)
    main(study_name='Serin_etal_2017_pheno', th=0)