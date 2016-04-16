

'''

INPUT: annotation file, ld file, result file
OUTPUT: Plot of posterior probabilities, plot of pvalues, ld plot, annotations color bar graph

'''

from string import letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv 
import matplotlib as mpl
from scipy.stats import norm
import math 

def Read_Data_CSV(filename):
	"""Reads in space delimited file as input and loads data into an array"""
	csv_file = csv.reader( open(filename, 'rb'),  delimiter = ' ' )
	#extract header line
	file_header = csv_file.next() 
	#read in all data from input file into n x 3 array
	data = [ row[:] for row in csv_file ] 
	data_array = np.array(data, dtype='double' )    
	np.set_printoptions(suppress=True) #suppress scientific notation  
	first_column = data_array[:,0]
	second_column = data_array[:,1]
	third_column = data_array[:,2]
	return data_array

def Read_Data_panda(filename):
	"""Reads in space delimted file as input and loads into a pandas data frame"""
	data = ld = pd.read_csv(filename, header=None, delimiter=r"\s+")
	return data

def Plot_Stats(filename):
	"""Reads in a space delimited result file and returns a scatterplot figure of the pvalues"""
	data = Read_Data_CSV(filename)
	position = data[:,0]
	zscore = data[:,1]
	#take -log10 scale of zscores 
	abs_zscore = np.absolute(zscore)
	pvalue = -1*(norm.logsf(abs_zscore)/math.log(10))
	plt.scatter(position, pvalue)
	plt.xlim(0,315)
	plt.ylim(0,52)
	plt.xlabel('SNPs')
	plt.ylabel('Observed association statistics (zscores)')
	plt.title('Plot of Zscores')
	plot_statistics = plt 
	return plot_statistics

def Plot_Probs(filename):
	"""Reads in a space delimited file and returns a scatterplot figure of the posterior probabilities"""
	data = Read_Data_CSV(filename)
	position = data[:,0]
	probabilities = data[:,2]
	plt.scatter(position, probabilities)
	plt.xlim(0,315)
	plt.ylim(0,.25)
	plt.xlabel('SNPs')
	plt.ylabel('Posterior Probabilities')
	plt.title('Plot of Posterior Probabilities')
	plot_probabilities = plt
	return plot_probabilities

def Plot_Correlation(filename):
	"""Reads in ld file and returns a correlation heatmap"""
	ld = Read_Data_panda(filename)
	sns.set(style="white")
	corr = ld.corr()
	mask = np.zeros_like(corr, dtype=np.bool)
	mask[np.triu_indices_from(mask)] = True
	#f, ax = plt.subplots(figsize=(11, 9))
	cmap = sns.diverging_palette(220, 10, as_cmap=True)
	ld_plot = sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3,
	           square=True, xticklabels=5, yticklabels=5,
	           linewidths=.5, cbar_kws={"shrink": .5})
	# Rotate and Crop figure
		#To be added 
	return ld_plot

def Plot_Annotations(filename):
	"""Reads in annotations file and returns a color bar graph"""
	csv_file = csv.reader( open(filename, 'rb'),  delimiter = ' ' )
	file_header = csv_file.next() 
	data = [ row[:] for row in csv_file ] 
	data_array = np.array(data, dtype='int_' ) 
	annotations = data_array[:,1]
	colors = []
	for a in annotations:
	    if a == 0:
	        colors.append('r')
	    else:
	        colors.append('b')
	fig = plt.figure(figsize=(8, 3))
	ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
	cmap = mpl.colors.ListedColormap(colors)
	cmap.set_over('0.25')
	cmap.set_under('0.75')
	N = len(annotations)
	bounds = range(1,N+1)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	annotations_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
	                                norm=norm,
	                                # to use 'extend', you must
	                                # specify two extra boundaries:
	                                #ticks=ticks,  #every 10, add a tick 
	                                spacing='proportional',
	                                orientation='horizontal')
	annotations_plot.set_label('Colorbar Test')
	annotations_plot = plt.figure()
	return annotations_plot


def main():
	filename = 'Enrich.DHS.38.annotations'
	plot = Plot_Annotations(filename)
	pyplot.show(plot)
	

if __name__ == "__main__": main()