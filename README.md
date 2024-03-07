# Motif Mark

## Overview
The goal of this algorithm is to visualize the location of DNA motifs on given DNA sequences to scale. It will account for motifs containing ambiguous nucleotide bases, and output a PNG image file that displays the location and size of each motif along the sequences given. 

## How to use the script
Before running the script, ensure your environment has ```pycairo``` installed. 

Run the python script with these inputs: 
* A fasta file containing gene sequences along which the motifs will be mapped. 
* A text file containing a list of motifs to be mapped.

The script uses argparse. Here is an example of how to run the script on the command line:

```./motif-mark-oop.py -f <path to fasta file> -m <path to motif txt file>```

## Notes
This algorithm takes a maximum of 5 different motifs because this is the maximum number of different colors that have been specified in the script. If your input contains more than 5 motifs, you will need to add additional RGB codes to the ```colors``` list in the script.  
