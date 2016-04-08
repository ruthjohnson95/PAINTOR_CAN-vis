import matplotlib.pyplot as plt
import numpy as np
import csv 
from matplotlib import pyplot
import matplotlib as mpl

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
fig = pyplot.figure(figsize=(8, 3))

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
plt.savefig('color_bar.png', bbox_inches='tight')
pyplot.show()



