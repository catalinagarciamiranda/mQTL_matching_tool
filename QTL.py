"""
This module provides functionalities for analyzing QTL (Quantitative Trait
Loci) data.

Author: Catalina Garcia Miranda

Date: 21/03/2024

Functions:
    get_region(markers, n):
        Retrieves the genomic region given a dataframe of markers and a marker index.

    Study:
        A class representing a QTL study.

        Methods:
            __init__(name):
                Initializes the Study instance with the given name and
                retrieves LOD scores and markers.

            get_peak(trait, threshold=0):
                Retrieves LOD scores and markers associated with a peak for a
                given trait.

            plot_scores(trait_id, threshold=0, marker=None, fig_title=None):
                Plots LOD scores for a given trait, optionally with a threshold
                line and markers.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams


def get_region(markers, n):
    # markers is df and n is a number tiki
    name = (markers.loc[n].name)
    chr = markers.loc[name].chr
    start_p = markers.loc[name].start - 50000
    end_p = start_p + 100000

    return chr, start_p, end_p

def find_local_maxima(data):    # NEW FUNCTION TO SEARCH LOCAL MAXIMUMS
    """Function to find local maxima in a dataset

    :param data: List or array containing the data
    :return: List of indices where local maxima occur
    """
    maxima_indices = []
    for i in range(1, len(data) - 1):
        if data[i] > data[i - 1] and data[i] > data[i + 1]:
            maxima_indices.append(i)
    return maxima_indices


class Study:

    def __init__(self, name):
        self.name = name

        base_url = 'https://www.bioinformatics.nl/AraQTLdev/media/data/' + \
                   self.name


        #if self.name == 'Joosen_etal_2013':
         #   path = f'{name}/{name}_LOD.txt'
          #  self.LOD_scores = pd.read_csv(path, sep="\t", index_col=0)
        #else:
            # retrieve LOD scores
        url_LOD = base_url + '/lod.txt'
        self.LOD_scores = pd.read_csv(url_LOD, sep="\t", index_col=0)


        # retrieve markers
        # retrieve marker
        url_marker = base_url + '/marker.txt'
        self.markers = pd.read_csv(url_marker, sep="\t", index_col=0)

    def get_peak(self, trait, threshold=0):
        """ Function to get the LOD and marker of a peak

        :param trait: the name of the metabolite of interest
        :param threshold: a number that determines the minimum score to be
        included as a peak

        :return: the LOD scores and markers associated with a peak
        """
        # key error:
        #print(trait in self.LOD_scores.index)
        if trait in self.LOD_scores.index:
            #print(trait)

            # get the LOD scores for a single trait
            trait_LOD = self.LOD_scores.loc[trait]

            # make all the values positive
            trait_LOD = abs(trait_LOD)

            # when no threshold is given, keep the 5% superior as peak
            if threshold == 0:
                threshold = trait_LOD.quantile(0.95)
                #print(threshold)

            peak_LOD = trait_LOD[trait_LOD >= threshold]
            labels = list(peak_LOD.index.astype(str))

            # keep the LOD scores above the threshold
            peak_markers = self.markers[self.markers.index.isin(labels)]
            if len(peak_markers.index) == 0:
                peak_markers = self.markers[self.markers.isin(labels).any(axis=1)]
            #print('peak_LOD = ', peak_LOD)
            #print('peak_markers = ', peak_markers)

            return peak_LOD, peak_markers

        else:
            return (None, None)


    def plot_scores(self, trait_id, threshold=0, marker=None, fig_title=None):

        # select the LOD scores for Arabidopsis trait_id
        trait_row = abs(self.LOD_scores.loc[trait_id])
        markers = self.markers

        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Gill Sans MT']

        # create the plot
        plt.figure(figsize=(12, 6))

        # plot LOD scores with a pastel blue color
        plt.plot(
            trait_row.index,
            trait_row.values,
            color='#25292B',
            label='LOD Scores'
        )

        if threshold == 0:
            threshold = trait_row.quantile(0.95)
            threshold = int(threshold*100)/100

        # plot the threshold line with a soft red color
        plt.axhline(
            threshold,
            color='#C67D08',
            linestyle='--',
            label=f'Threshold {threshold}',
            alpha=0.9
        )

        if marker is not None:
            colormap = plt.cm.viridis    # set the color palette for the lines
            idx = 1
            for m in marker:
                color = colormap(idx/len(marker))
                idx += 1
                # Add a vertical line or marker for the specified marker index
                plt.axvline(
                    m,
                    linestyle='--',
                    label=f'Marker {m}',
                    color=color,
                    alpha=0.8   # transparency
                )


        if fig_title == None:
            fig_title = trait_id

        plt.title(
            f'QTL profile for trait: {fig_title}, '
            f'threshold = {threshold}'
        )
        plt.xlabel('Markers')

        # Show every other x-axis label
        increment = 20
        plt.xticks(
            np.arange(len(trait_row.index))[::increment],
            trait_row.index[::increment],
            rotation=45,
            ha='right'
        )

        plt.ylabel('LOD score')

        # alternate background color to indicate chromosomes
        for i, m in enumerate(markers.index):

            if markers.loc[m].chr % 2 == 0:
                # light gray background for even chromosomes
                color = '#D8D8D8'
            else:
                # light beige background for odd chromosomes
                color = '#F0EAD6'

            plt.axvspan(
                i - 0.5,
                i + 0.5,
                facecolor=color,
                alpha=0.5
            )

        plt.legend()
        plt.tight_layout()
        #plt.show(block=False)

    def get_peak_local_max(self, trait, threshold=0):
        """ Function to get the LOD and marker of local maxima peaks

        :param trait: the name of the metabolite of interest
        :param threshold: a number that determines the minimum score to be
        included as a peak

        :return: the LOD scores and markers associated with local maxima peaks
        """
        if trait in self.LOD_scores.index:
            # Get the LOD scores for a single trait
            trait_LOD = self.LOD_scores.loc[trait]

            # Make all the values positive
            trait_LOD = abs(trait_LOD)

            # When no threshold is given, keep the 5% superior as peak
            if threshold == 0:
                threshold = trait_LOD.quantile(0.95)

            # Find local maxima
            maxima_indices = find_local_maxima(trait_LOD)

            # Filter maxima by threshold
            peak_indices = [i for i in maxima_indices if trait_LOD[i] >= threshold]

            # Get the LOD scores and markers associated with the peak indices
            peak_LOD = trait_LOD.iloc[peak_indices]
            peak_markers = self.markers.iloc[peak_indices]

            return peak_LOD, peak_markers
        else:
            return None, None


