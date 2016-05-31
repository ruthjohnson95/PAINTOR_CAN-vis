# INSTRUCTIONS

1. Install additional libraries by running `install_mac.sh` or `install_linux.sh`.
2. Run the following command inside of the CANVIS_sample folder 
```
$ python CANVAS.py -l chr4.3473139.rs6831256.pleio.annot -z ldl.Zscore tg.Zscore -a chr4.3473139.rs6831256.annotations -s E066-H3K27ac.narrowPeak.Adult_Liver E066-H3K4me1.narrowPeak.Adult_Liver -r chr4.3473139.rs6831256.ld  -t 99
```
The image produced should match the `example_fig.pdf` file given in the sample. 