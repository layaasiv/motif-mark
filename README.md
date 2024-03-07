# Motif Mark

## Overview
The goal of this algorithm is to visualize the location of DNA motifs on given DNA sequences to scale. It will account for motifs containing ambiguous nucleotide bases, and output a PNG image file that displays the location and size of each motif along the sequences given. 

## Usage
Before running the script, ensure your environment has ```pycairo``` installed. 

Run the python script with these inputs: 
* A fasta file containing gene sequences along which the motifs will be mapped. 
* A text file containing a list of motifs to be mapped.

The script uses argparse with 2 required arguments: ```-f``` and ```-m```. Here is an example of how to run the script on the command line:

```./motif-mark-oop.py -f <path to fasta file> -m <path to motif txt file>```

The algorithm only takes 1 fasta and motif file at a time, and will depict motif locations for each sequence in that file in the image. The final PNG file corresponding to the inputted fasta file will have the same filename as the fasta file.

## Notes
This algorithm takes a maximum of 5 different motifs because this is the maximum number of different colors that have been specified in the script. If the ```-m``` input contains more than 5 motifs, additional RGB codes must be added to the ```colors``` list in the script. 

Lines represent introns, and black rectangles represent exons. Rectangles of other colors represent different motifs (specified in the legend of the image).

## Packages used 
* argparse
* re
* itertools
* cairo
