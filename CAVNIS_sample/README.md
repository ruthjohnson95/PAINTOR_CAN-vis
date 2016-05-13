# INSTRUCTIONS

1. Install additional libraries by running `install.sh`
2. Run the following command inside of the CANVIS_sample folder 
```
python CANVIS_beta.py -l chr4.3473139.rs6831256.post.filt.300 -z ldl.Zscore -r chr4.3473139.rs6831256.ld.filt.300 -a chr4.3473139.rs6831256.annot.filt.300 -s E066.H3K27ac.narrowPeak.Adult_Liver E066.H3K4me1.narrowPeak.Adult_Liver -t 99 -g y
```
3. A svg file named fig_final.svg should be produced and match the example_fig.pdf file