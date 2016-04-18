from string import letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import matplotlib as mpl
from scipy.stats import norm
import math
from optparse import OptionParser
import svgutils.transform as sg
import sys

def Read_Input(locus_fname, ld_fname, annotation_fname):
    """Function that reads in all your data """
    csv_file = csv.reader(open(locus_fname, 'rb'), delimiter=' ')
    file_header = csv_file.next()  # extract header line
    locus_data = [row[:] for row in csv_file]
    locus = np.array(locus_data, dtype='double')
    ld = pd.read_csv(ld_fname, header=None, delimiter=r"\s+")
    csv_file = csv.reader(open(annotation_fname, 'rb'), delimiter=' ')
    file_header = csv_file.next()
    annotation_data = [row[:] for row in csv_file]
    annotation = np.array(annotation_data, dtype='int_')
    annotation = annotation[:, 1]
    return [locus, ld, annotation]

def Plot_Position_Value(position, zscore, pos_prob ):

    """Function that plots z-scores, posterior probabilites, other features """
    fig = plt.figure(figsize=(12, 6.25))
    sub1 = fig.add_subplot(2, 1, 1, axisbg='white')
    pvalue = Zscore_to_Pvalue(zscore)
    sub1.scatter(position, pvalue)
    sub2 = fig.add_subplot(2, 1, 2, axisbg='white')
    sub2.scatter(position, pos_prob)
    value_plots = fig
    return value_plots #returns subplots with both graphs

def Plot_Heatmap(correlation_matrix):
    """Function that plots heatmap of LD matrix"""
    fig = plt.figure(figsize=(12, 6.25))
    sns.set(style="white")
    correlation = correlation_matrix.corr()
    mask = np.zeros_like(correlation, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(correlation, mask=mask, cmap=cmap, vmax=.3 ,square=True, xticklabels=5,
                yticklabels=5, linewidths=.5, cbar_kws ={"shrink": .5})
    heatmap = fig
    # Rotate and Crop figure to be added later
    return heatmap

def Plot_Annotations(annotation_vectors, annotation_names):
    """Plot the annotations with labels"""
    colors = []
    for a in annotation_vectors:
        if a == 0:
            colors.append('r')
        else:
            colors.append('b')
    fig = plt.figure(figsize=(12, 2))
    ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
    cmap = mpl.colors.ListedColormap(colors)
    cmap.set_over('0.25')
    cmap.set_under('0.75')
    N = len(annotation_vectors)
    bounds = range(1, N+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    annotation_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                    norm=norm,
                                    spacing='proportional',
                                    orientation='horizontal')
    annotation_plot.set_label('Annotations')
    annotation_plot = plt
    return annotation_plot

def Assemble_Figure(value_plots, heatmap, annotation_plot):
    """Assemble everything together"""
    value_plots.savefig('value_plots.svg', format='svg', dpi=1200)
    heatmap.savefig('heatmap.svg', format='svg', dpi=1200)
    annotation_plot.savefig('annotation_plot.svg', format='svg', dpi=1200)

    fig = sg.SVGFigure("13in", "19in")
    value_plots = sg.fromfile('value_plots.svg')
    heatmap  = sg.fromfile('heatmap.svg')
    annotation_plot = sg.fromfile('annotation_plot.svg')

    plot1 = value_plots.getroot()
    plot2 = heatmap.getroot()
    plot3 = annotation_plot.getroot()
    plot2.moveto(0, 450, scale=1.2)
    plot3.moveto(0, 950, scale=.8)

    fig.append([plot1, plot2, plot3])
    fig.save("fig_final.svg")


def Zscore_to_Pvalue(zscore):
    """Function that converts zscores to pvalues"""
    abs_zscore = np.absolute(zscore)
    pvalue = -1 * (norm.logsf(abs_zscore) / math.log(10))
    return pvalue

def main():

    # defaults
    plot_annotations = None

    # Parse the command line data
    parser = OptionParser()
    parser.add_option("-l", "--locus_name", dest="locus_name")
    parser.add_option("-a", "--annotation_name", dest="annotation_name")
    parser.add_option("-r", "--ld_name", dest="ld_name")
    parser.add_option("-p", "--plot_annotations", dest="plot_annotations")
    # add other optinos to parse

    # extract options
    (options, args) = parser.parse_args()

    locus_name = options.locus_name
    annotation_name = options.annotation_name
    ld_name = options.ld_name
    plot_annotations = options.plot_annotations
    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position) *
        --ld_name [r] specify the ld_matrix file name *
        --annotation_name [-a]  specify annotation file name *
        --plot_annotations [-p] specify which annotations to plot [default: None]
        """

    # check if required flags are presnt
    """if(locus_name == None or annotation_name == None or ld_name == None or plot_annotations == None):
        sys.exit(usage)"""

    [locus, ld, annotation] = Read_Input(locus_name, ld_name, annotation_name)
    value_plots = Plot_Position_Value(locus[:, 0], locus[:, 1],locus[:, 2])
    heatmap = Plot_Heatmap(ld)
    annotation_plot = Plot_Annotations(annotation, annotation_name)

    Assemble_Figure(value_plots, heatmap, annotation_plot)



if __name__ == "__main__":
    main()
