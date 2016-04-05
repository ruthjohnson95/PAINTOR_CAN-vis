#data.py

import csv 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt
import math 
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
plt.scatter(pos_array, pvalue)
plt.xlim(0,315)
plt.ylim(0,52)
plt.xlabel('SNPs')
plt.ylabel('Observed association statistics (zscores)')
plt.title('Plot of Zscores')
plt.show()

#scatterplot of probabilities 
plt.scatter(pos_array, prob_array)
plt.xlim(0,315)
plt.ylim(0,.25)
plt.xlabel('SNPs')
plt.ylabel('Posterior Probabilities')
plt.title('Plot of Posterior Probabilities')
plt.show()




