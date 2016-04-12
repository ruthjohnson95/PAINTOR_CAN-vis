"""Examples illustrating the use of plt.subplots().

This function creates a figure and a grid of subplots with a single call, while
providing reasonable control over how the individual plots are created.  For
very refined tuning of subplot creation, you can still use add_subplot()
directly on a new figure.
"""
import csv 
import numpy as np 
import math 
import seaborn as sns 
import matplotlib.pyplot as plt
#from matplotlib import pyplot
import matplotlib as mpl
from scipy.stats import norm

filename = 'Enrich_DHS_4_result'
#read input file 
csv_file = csv.reader( open(filename, 'rb'),  delimiter = ' ' )
#extract header line
file_header = csv_file.next() 

#read in all data from input file into n x 3 array
data = [ row[:] for row in csv_file ] #gets last row 
data_array = np.array(data, dtype='double' )    
#print(data_array) #testing purposes 
np.set_printoptions(suppress=True) #suppress scientific notation 
#position array
print "Printing positions"
pos_array = data_array[:,0]
print(pos_array)
#zscore array
print "Printing zscores"
zscore_array = data_array[:,1]
print(zscore_array)
#posterior probability array
print "Printing posterior probabities"
prob_array = data_array[:,2]
print(prob_array)

#take -log10 
abs_zscore_array = np.absolute(zscore_array)
pvalue = -1*(norm.logsf(abs_zscore_array)/math.log(10)) #convert to base 10

#plot scatterplot 

#plt.show()


# Simple data to display in various forms
x = np.linspace(0, 100,100)
y = np.sin(x ** 2)

plt.close('all')

#read in file 
filename = 'Enrich.DHS.38.annotations'
csv_file = csv.reader( open(filename, 'rb'),  delimiter = ' ' )
#extract header line
file_header = csv_file.next() 
data = [ row[:] for row in csv_file ] #gets last row 
data_array = np.array(data, dtype='int_' )  

annotations = data_array[:,1]


#annotations = [0,0,1,0,1,1,1,0,1]
colors = []

for a in annotations:
    if a == 0:
        colors.append('r')
    else:
        colors.append('b')  

print colors 

# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(8, 3))

ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])

# The second example illustrates the use of a ListedColormap, a
# BoundaryNorm, and extended ends to show the "over" and "under"
# value colors.
cmap = mpl.colors.ListedColormap(colors)
cmap.set_over('0.25')
cmap.set_under('0.75')

# If a ListedColormap is used, the length of the bounds array must be
# one greater than the length of the color list.  The bounds must be
# monotonically increasing.
N = len(annotations)
bounds = range(1,N+1) #note, outer bound must be n+1 
print bounds 

#ticks = np.arange(1,N+1,10)
#print ticks

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                norm=norm,
                                # to use 'extend', you must
                                # specify two extra boundaries:
                                #ticks=ticks,  #every 10, add a tick 
                                spacing='proportional',
                                orientation='horizontal')
cb2.set_label('Colorbar Test')
#plt.savefig('color_bar.png', bbox_inches='tight')
#pyplot.show()


# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(3, sharex=True)
axarr[0].scatter(pos_array, prob_array)
axarr[0].set_title('Sharing X axis')
axarr[1].scatter(pos_array, pvalue)

#axarr[2].cb2
'''
plt.xlim(0,315)
plt.ylim(0,52)
plt.xlabel('SNPs')
plt.ylabel('Observed association statistics (zscores)')
plt.title('Plot of Zscores')
'''

plt.show()