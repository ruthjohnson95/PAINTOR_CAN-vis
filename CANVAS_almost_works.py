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
import cairosvg
import sys


def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []

    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

def Read_Input(locus_fname, ld_fname, annotation_fname):
    """Function that reads in all your data """
    csv_file = csv.reader(open(locus_fname, 'rb'), delimiter=' ')
    file_header = csv_file.next()  # extract header line
    locus_data = [row[:] for row in csv_file]
    locus = np.array(locus_data, dtype='double')
    ld = pd.read_csv(ld_fname, header=None, delimiter=r"\s+")
    annotation = []
    for i in range(0, len(annotation_fname)):
        csv_file = csv.reader(open(annotation_fname[i], 'rb'), delimiter=' ')
        file_header = csv_file.next()
        annotation_data = [row[:] for row in csv_file]
        annotation_array = np.array(annotation_data, dtype='int_')
        annotation_array = annotation_array[:, 1]
        annotation.append(annotation_array)
    return [locus, ld, annotation]

def Plot_Position_Value(position, zscore, pos_prob ):

    """Function that plots z-scores, posterior probabilites, other features """
    fig = plt.figure(figsize=(12, 6.25))
    sub1 = fig.add_subplot(2, 1, 1, axisbg='white')
    plt.xlim(np.amin(position), np.amax(position))
    plt.ylabel('-log10(pvalue)')
    pvalue = Zscore_to_Pvalue(zscore)
    sub1.scatter(position, pvalue)
    sub2 = fig.add_subplot(2, 1, 2, axisbg='white')
    plt.xlim(np.amin(position), np.amax(position))
    plt.gca().set_ylim(bottom=0)
    plt.ylabel('Posterior probabilities')
    plt.xlabel('Location')
    sub2.scatter(position, pos_prob)
    value_plots = fig
    return value_plots #returns subplots with both graphs

def Plot_Heatmap(correlation_matrix, hue1, hue2):
    """Function that plots heatmap of LD matrix"""
    fig = plt.figure(figsize=(6.25, 6.25))
    sns.set(style="white")
    correlation = correlation_matrix.corr()
    mask = np.zeros_like(correlation, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    h1 = int(hue1)
    h2 = int(hue2)
    cmap = sns.diverging_palette(h1, h2, as_cmap=True)
    sns.heatmap(correlation, mask=mask, cmap=cmap, vmax=.3, square=True,
                linewidths=.5, cbar=False, xticklabels=False, yticklabels=False, ax=None)
    heatmap = fig
    return heatmap

def Plot_Annotations(annotation_names, annotation_vectors):
    """Plot the annotations with labels"""
    annotation_tuple = []
    for i in range(0, len(annotation_names)):
        annotation = annotation_vectors[i]
        colors = []
        for a in annotation:
            if a == 0:
                colors.append('#DADFE1')
            else:
                colors.append('#4B77BE')
        fig = plt.figure(figsize=(12, 1.0))
        ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
        cmap = mpl.colors.ListedColormap(colors)
        cmap.set_over('0.25')
        cmap.set_under('0.75')
        n = len(annotation)
        bounds = range(1, n+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        annotation_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                                     norm=norm,
                                                     spacing='proportional',
                                                     orientation='horizontal')
        annotation_plot.set_label(annotation_names[i])
        annotation_plot = plt
        annotation_tuple.append(annotation_plot)
    return annotation_tuple

def Assemble_Figure(value_plots, heatmap, annotation_plot):
    """Assemble everything together"""
   # value_plots.savefig('value_plots.svg', format='svg', dpi=1200)
   # heatmap.savefig('heatmap.svg', format='svg', dpi=1200)
    fig = sg.SVGFigure("13in", "23in")
   # value_plots = sg.fromfile('value_plots.svg')
   # heatmap = sg.fromfile('heatmap.svg')
   # plot1 = value_plots.getroot()
   # plot2 = heatmap.getroot()
   # plot2.rotate(-45, 600, 600)
   # plot2.moveto(250, 150, scale=1.375)
   # fig.append([plot2, plot1])
    #annotation plots
    for i in annotation_plot:
        plot = i
        plot.savefig('annotation_plot.svg', format='svg', dpi=1200)
        plot3 = annotation_plot.getroot()
        plot3.move(60, 300 + 72*(i+1), scale=.9) #CAN'T USE I HERE B/C IT'S AN OBJECT; ENUMERATE INSTEAD
        fig.append(plot3)

    fig.save("fig_final.svg")
    cairosvg.svg2pdf(url='fig_final.svg', write_to='fig_final.pdf')


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
    #parser.add_option("-a", "--annotation_name", dest="annotation_name")
    parser.add_option("-n", "--number_of_args", )
    parser.add_option("-a", "--annotation_names", dest="annotation_names", action='callback', callback=vararg_callback)
    parser.add_option("-p", "--annotation_plot", dest="annotation_plot", action='callback', callback=vararg_callback)
    parser.add_option("-r", "--ld_name", dest="ld_name")
    parser.add_option("--h1", "--hue1", dest="hue1", default=233)
    parser.add_option("--h2", "--hue2", dest="hue2", default=155)

    # extract options
    (options, args) = parser.parse_args()
    locus_name = options.locus_name
    ld_name = options.ld_name
    annotation_names = options.annotation_names
    annotation_plot = options.annotation_plot
    hue1 = options.hue1
    hue2 = options.hue2
    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position) *
        --ld_name [r] specify the ld_matrix file name *
        --annotation_name [-a]  specify annotation file name *
        --plot_annotations [-p] specify which annotations to plot [default: None]
        --hue1 [-h1]
        --hue2 [-h2]
        """

    # check if required flags are presnt
    """if(locus_name == None or annotation_name == None or ld_name == None or plot_annotations == None):
        sys.exit(usage)"""

    [locus, ld, annotation] = Read_Input(locus_name, ld_name, annotation_plot)
    value_plots = Plot_Position_Value(locus[:, 0], locus[:, 1], locus[:, 2])
    heatmap = Plot_Heatmap(ld, hue1, hue2)
    annotation_plot = Plot_Annotations(annotation_names, annotation)

    Assemble_Figure(value_plots, heatmap, annotation_plot)



if __name__ == "__main__":
    main()
