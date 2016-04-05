from string import letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image

sns.set(style="white")

#read in input files 
filename = 'Enrich.DHS.6.ld'
ld = pd.read_csv(filename, header=None, delimiter=r"\s+")

#create correlation matrix
corr = ld.corr()

# Generate a mask for the upper triangle
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3,
            square=True, xticklabels=5, yticklabels=5,
            linewidths=.5, cbar_kws={"shrink": .5}, ax=ax)

plt.savefig('fig.png', bbox_inches='tight')
im1 = Image.open("fig.png")
im2 = im1.rotate(45, expand=True)
im2 = im2.crop((0, 0, 10, 10)
#im2.show()
im2.save("fig_rotate.png")
#plt.show()
