# CANVAS (Correlation ANnotation Visualization Association Statistics)

## Introduction
We provide here documentation for CANVAS, a visualization tool utilized for the output of PAINTOR. We will describe the  neccessary libraries, file formats, and input parameters, as well as a sample script. 

## File Formats

#### Zscore File
The zscore file is space delimted with the first row is expected to be a header row containing the labels of each column. This file must contain at least one column of zscores with the name of the specific zscore as the heading, one column of posterior probabilites with the heading `Posterior_Prob`, and a column holding the positions with the heading `pos`. CANVAS can plot between 1 - 3 zscores. An example of a correctly formatted file is shown below: 

```
chr pos rsid hdl.A0 hdl.A1 hdl.Zscore ldl.Zscore tc.Zscore tg.Zscore Posterior_Prob
chr4 3349626 rs74823150 A G -0.594734 0.696424 0.507638 0.345514 2.7728e-06
chr4 3350248 rs2749779 A G -1.320755 1.068966 1.754386 2.153846 1.87298e-06

```
#### Annotation File
The annotation file is space delimited, and should contain a column per annotation where each position has a 0 or 1 to represent a specific annotation. The first row is assumed to be a header row containing the names of the annotations. The file contain multiple columns, but a max of 5 annotations can be plotted. 

```
E066-H3K27ac.narrowPeak.Adult_Liver E066-H3K4me1.narrowPeak.Adult_Liver
0 1
0 0
1 0
```
#### LD file
The ld matrix file is a space delimited file with no header row or indexing column. Each row should contain the correlations corresponding to other locations of the sample. 

```
0.20794 -0.33251 -0.0018495 -0.33251 0.7164 -0.33181 0.70725 0.62996 -0.33181 -0.33181 -0.33181
```

## Neccessry Installed Libraries 
```
string
numpy
panda
seaborn
matplotlib
scipy
math
optparse
svgutils
cairosvg
sys

```
svgutils is not included in the Anaconda distribution, and can be installed with following: `pip install svgutils --user`

cariosvg is also not included in the Anaconda distribution, and can be installed with following: `brew install cairo` and then `$ pip install cairosvg`

## Input 
```
--locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position) 
--ld_name [r] specify the ld_matrix file name
--annotation_name [-a]  specify annotation file name
--plot_annotations [-p] specify which annotations to plot [default: None]
--hue1 [-h1] (degreees, 0 - 360)
--hue2 [-h2] (degrees, 0 - 360)

```
When including 

### Output
1. Scatterplot of location versus -log10(pvalue)
2. Scatterplot of location versus posterior probabilites
3. Correlation heatmap 
4. Annotation bars with each annotation defined as a different color 

### Example Script 
```
python canvas.py -l [path to locus file] -r [path to ld matrix file] -a [path to annottions file] -p [annotation names] -h1 [0-360] -h2 [0-360]

```



