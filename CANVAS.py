import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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

def Read_Input(locus_fname, zscore_names, ld_fname, annotation_fname, specific_annotations):
    """Function that reads in all your data """
    zscore_data = pd.read_csv(locus_fname, delim_whitespace=True)
    zscores = zscore_data[zscore_names]
    location = zscore_data['pos']
    pos_prob = zscore_data['Posterior_Prob']
    ld = pd.read_csv(ld_fname, header=None, delim_whitespace=True)
    annotation_data = pd.read_csv(annotation_fname, delim_whitespace=True)
    annotations = annotation_data[specific_annotations]
    zscores = zscores.as_matrix()
    pos_prob = pos_prob.as_matrix()
    location = location.as_matrix()
    annotations = annotations.as_matrix()
    return [zscores,pos_prob,location, ld, annotations]

def Plot_Statistic_Value(position, zscore, zscore_names):
    zscore_tuple = []
    color_array = ['#D64541', '#1E824C', '#F89406']
    for i in range(0, len(zscore_names)):
        fig = plt.figure(figsize=(12, 3.75))
        sub = fig.add_subplot(1,1,1, axisbg='white')
        plt.xlim(np.amin(position), np.amax(position) + 1)
        plt.ylabel('-log10(pvalue)')
        plt.xlabel(zscore_names[i])
        z = zscore[:, i]
        pvalue = Zscore_to_Pvalue(z)
        sub.scatter(position, pvalue, color=color_array[i])
        plt.gca().set_ylim(bottom=0)
        value_plot = fig
        zscore_tuple.append(value_plot)
    return zscore_tuple

def Plot_Position_Value(position, pos_prob):
    """Function that plots z-scores, posterior probabilites, other features """
    [credible_loc, credible_prob] = Credible_Set(position, pos_prob, .99)
    fig = plt.figure(figsize=(12, 3.25))
    sub1 = fig.add_subplot(1,1,1, axisbg='white')
    plt.xlim(np.amin(position), np.amax(position)+1)
    plt.ylabel('Posterior probabilities')
    plt.xlabel('Location')
    sub1.scatter(position, pos_prob, color='#2980b9', label='Non-Credible Set')
    sub1.scatter(credible_loc, credible_prob, color='#D91E18', label='Credible Set')
    legend = plt.legend(loc='upper right', shadow=True)
    frame = legend.get_frame()
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    plt.gca().set_ylim(bottom=0)
    value_plots = fig
    return value_plots

def Credible_Set(position, pos_prob, threshold):
    total = sum(pos_prob)
    bounds = threshold*total
    #make into tuples
    tuple_vec = []
    for i in range(0, len(position)):
        tup = (position[i], pos_prob[i])
        tuple_vec.append(tup)
    #order tuple from largest to smallest
    tuple_vec = sorted(tuple_vec, key=lambda x: x[1], reverse=True)
    credible_set_value = []
    credible_set_loc = []
    total = 0
    for tup in tuple_vec:
        total += tup[1]
        credible_set_loc.append(tup[0])
        credible_set_value.append(tup[1])
        if total > bounds:
            break
    return credible_set_loc, credible_set_value

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
    sns.heatmap(correlation, mask=mask, cmap=cmap, square=True,
                linewidths=0, cbar=False, xticklabels=False, yticklabels=False, ax=None)
    heatmap = fig
    return heatmap

def Plot_Annotations(annotation_names, annotation_vectors):
    """Plot the annotations with labels"""
    annotation_tuple = []
    color_array = ['#663399', '#e74c3c', '#049372', '#F89406', '#1E8BC3']
    for i in range(0, len(annotation_names)):
        annotation = annotation_vectors[:,i]
        colors = []
        for a in annotation:
            if a == 1:
                colors.append(color_array[i])
            else:
                colors.append('#ecf0f1')
        fig = plt.figure(figsize=(12, 1.0))
        ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
        cmap = mpl.colors.ListedColormap(colors)
        cmap.set_over('0.25')
        cmap.set_under('0.75')
        bounds = range(1, len(annotation)+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        annotation_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                                     norm=norm,
                                                     spacing='proportional',
                                                     orientation='horizontal')
        annotation_plot.set_label(annotation_names[i])
        annotation_plot.set_ticks([])
        annotation_plot = fig
        annotation_tuple.append(annotation_plot)
    return annotation_tuple


def Assemble_Figure(stats_plot, value_plots, heatmap, annotation_plot):
    """Assemble everything together"""
    fig = sg.SVGFigure("9in", "14in")
    value_plots.savefig('value_plots.svg', format='svg', dpi=1200)
    heatmap.savefig('heatmap.svg', format='svg', dpi=1200)
    value_plots = sg.fromfile('value_plots.svg')
    heatmap = sg.fromfile('heatmap.svg')
    plot1 = value_plots.getroot()
    plot4 = heatmap.getroot()

    #transform and add heatmap figure
    y_scale = 50 * (len(annotation_plot)) + len(stats_plot)*275 + 275
    if len(annotation_plot) == 1:
        y_scale = 0
    plot4.moveto(-10, y_scale, scale=1.425)
    plot4.rotate(-45, 0, 0)
    fig.append(plot4)

    #transform and add value plot
    y_move = 275*len(stats_plot)
    plot1.moveto(0, y_move)
    fig.append(plot1)

    #transform and add zscore plots
    index = 0
    for plot in stats_plot:
        plot.savefig('stats_plot.svg', format='svg', dpi=1200)
        plot = sg.fromfile('stats_plot.svg')
        plot2 = plot.getroot()
        y_move = 275 * index
        index += 1
        plot2.moveto(0, y_move)
        fig.append(plot2)

    #transform and add annotations plots
    index = 0
    for plot in annotation_plot:
        plot.savefig('annotation_plot.svg', format='svg', dpi=1200)
        plot = sg.fromfile('annotation_plot.svg')
        plot3 = plot.getroot()
        y_move = 200+50*(index+1) + 275*len(stats_plot)
        plot3.moveto(60, y_move, scale=.9)
        index += 1
        fig.append(plot3)
    #export final figure as a svg and pdf
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
    parser.add_option("-z", "--zscores", dest="zscores", action='callback', callback=vararg_callback)
    parser.add_option("-a", "--annotations", dest="annotations")
    parser.add_option("-s", "--specific_annotations", dest="specific_annotations", action='callback', callback=vararg_callback)
    parser.add_option("-r", "--ld_name", dest="ld_name")
    parser.add_option("--h1", "--hue1", dest="hue1", default=240)
    parser.add_option("--h2", "--hue2", dest="hue2", default=10)

    # extract options
    (options, args) = parser.parse_args()
    locus_name = options.locus_name
    zscore_names = options.zscores
    ld_name = options.ld_name
    annotations = options.annotations
    annotation_names = options.specific_annotations
    hue1 = options.hue1
    hue2 = options.hue2
    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position) *
        --ld_name [r] specify the ld_matrix file name *
        --annotations [-a]  specify annotation file name *
        --s [-s] specify which annotations to plot [default: None]
        --hue1 [-h1]
        --hue2 [-h2]
        """

    #check if required flags are presnt
    #if(locus_name == None or annotation_names == None or ld_name == None):
        #sys.exit(usage)

    [zscores, pos_prob, location, ld, annotations] = Read_Input(locus_name, zscore_names, ld_name, annotations, annotation_names)
    stats_plot = Plot_Statistic_Value(location, zscores, zscore_names)
    value_plots = Plot_Position_Value(location, pos_prob)
    heatmap = Plot_Heatmap(ld, hue1, hue2)
    annotation_plot = Plot_Annotations(annotation_names, annotations)

    Assemble_Figure(stats_plot, value_plots, heatmap, annotation_plot)



if __name__ == "__main__":
    main()
