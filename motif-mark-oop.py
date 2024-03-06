#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
import cairo
import itertools

# set up argpare
def get_args():
    parser = argparse.ArgumentParser(description="a")
    parser.add_argument("-f", "--fastafile", help="file path of input FASTA file containing sequences to be parsed for presence of motifs.", required=True)
    parser.add_argument("-m", "--motiffile", help="file path of input text file containing all motifs to be searched for.", required=True)
    return parser.parse_args()

args = get_args()
fasta_file = args.fastafile
motifs_file = args.motiffile

# define classes
class Motif:
    def __init__(self, the_seq: str, the_color: tuple, the_len: int, legend_y: int):
        '''Definition of a Motif class.'''
        ## Data ##
        self.seq = the_seq
        self.color = the_color
        self.length = the_len
        self.legend_ypos = legend_y
    ## Methods ##
    def all_possible_seq(self, seq: str) -> list:
        '''Gives a list of all possible motif sequences for those containing ambiguous bases.'''
        ambig = {"Y": ["T", "C"], 
                 "V":["A", "C", "G"], 
                 "R":["G", "A"], 
                 "M":["A", "C"],
                 "K":["G", "T"], 
                 "S":["G", "C"], 
                 "W":["A", "T"], 
                 "H":["A", "C", "T"], 
                 "B":["G", "T", "C"],
                 "D":["G", "A", "T"], 
                 "N":["A", "C", "T", "G"], 
                 "U":["T"]
                }
        seq = seq.upper()
        groups = itertools.groupby(seq, lambda char:char not in ambig)

        bases_per_pos = []
        for key, group in groups:
            if key:
                bases_per_pos.extend([[nuc] for nuc in group])
            else:
                for nuc in group:
                    bases_per_pos.append(ambig[nuc])
        all_seq = [''.join(base) for base in itertools.product(*bases_per_pos)]
        return all_seq

    def motif_finder(self, gene1_seq: str, the_seq: str) -> list:
        '''Finds start and end index position of a given motif in a given gene sequence. The search is case
        insensitive.'''
        match_pos = []
        reg = re.finditer(the_seq, gene1_seq, re.IGNORECASE)
        for match in reg:
            pos = match.span()
            match_pos.append(pos)
        return match_pos
    
    def all_motif_positions(self, motif1_all: list, gene1_seq: str) -> list:
        '''Searches the gene sequence for all possible instances of a motif (motif1_all is a list of all possible
        sequences the motif could present as if it has ambiguous bases). Returns a list of lists, each containing
        tuples of the start and end index positions of a motif.'''
        motif_pos = []
        for motif in motif1_all:
            motif_pres = motif1.motif_finder(gene1_seq, motif)
            if motif_pres != []:
                motif_pos.append(motif_pres)
            else:
                pass
        return motif_pos

    def unlisting_fxn(self, match_pos: list) -> list:
        '''Takes the product of the all_motif_positions function and combines all tuples into the same list. Returns
        list of tuples of the start and end index positions of a motif.'''
        motif_pos_unlisted = []
        if match_pos != []:
            for lst in match_pos:
                for coordinate in lst:
                    motif_pos_unlisted.append(coordinate)
        return motif_pos_unlisted 

class Gene:
    def __init__(self, the_name: str, the_len: int, the_seq: str, the_diag_yposition: int, the_name_yposition: int):
        '''Definition of a Gene class.'''
        ## Data ##
        self.name = the_name
        self.seq = the_seq
        self.len = the_len
        self.diag_ypos = the_diag_yposition
        self.name_ypos = the_name_yposition
    ## Methods ##
    def exon_finder(self, the_seq: str) -> list:
        '''Finds start and end index positions of exons in a gene sequence given as a string, assuming exons are 
        denoted as capitalized letters and introns are lowercase letters.'''
        match_pos = []
        matches  = re.finditer("[A-Z]+", the_seq)
        for match in matches:
            pos = match.span()
            match_pos.append(pos)
        return match_pos

class Draw:
    def __init__(self, init_xpos: int, init_ypos: int, canvas_width: int, canvas_height: int):
        '''Definition of a Draw class.'''
        ## Data ##
        self.init_xpos = init_xpos
        self.init_ypos = init_ypos
        self.canv_wid = canvas_width
        self.canv_hgt = canvas_height

    ##Methods##
    def init_canvas(self, init_xpos: int, init_ypos: int, canv_wid: int, canv_hgt: int):
        '''This function initializes the canvas and context. Input canvas width and height. Also input the initial
        position of pointer on the canvas.'''
        surface = cairo.PDFSurface(output_png + ".pdf", canv_wid, canv_hgt)
        context = cairo.Context(surface)
        context.save()
        context.set_source_rgb(1, 1, 1)
        context.paint()
        context.restore()
        context.move_to(init_xpos,init_ypos)
        return surface, context
    
    def write_gene_name(self, context, gene_name: str, xpos: int, ypos: int):
        '''This function writes the gene name on the canvas. Input gene name and x,y position of where the name should
        be written.'''
        context.set_source_rgb(0,0,0)
        context.move_to(xpos, ypos) 
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.set_font_size(16)
        context.show_text(gene_name)
        return context

    def draw_gene(self, context, xpos: int, gene_ypos: int, gene_len: int):
        '''This function draws a line (width 4) to represent the gene length. Input x,y position where gene should
        be drawn. Also input the base length of the gene sequence.'''
        context.set_source_rgb(0,0,0)
        context.set_line_width(4)
        context.move_to(xpos,gene_ypos)
        context.line_to(gene_len + 100, gene_ypos)
        context.stroke()
        return context
    
    def draw_exon(self, context, exon_positions: list, gene_ypos: int):
        '''This function draws the exons of the gene sequence as black rectangles. Input the exon positions and y 
        poition corresponding to the gene diagram location on the canvas.'''
        for exon_coordinate in exon_positions:
            exon_len = exon_coordinate[1] - exon_coordinate[0]

            context.rectangle(exon_coordinate[0] + 101,
                            gene_ypos - 10,
                            exon_len,
                            20)
            context.fill()
        return context
    
    def draw_motif(self, context, motif_positions: list, motif_color: tuple, gene_ypos: int, motif_len: int):
        '''This function draws the motifs on the gene sequence. Input the motif start/end positions as a list, RGB 
        value of motif color as a tuple. Also input y position corresponding to the gene diagram location, and motif
        base pair length.'''
        if motif_positions != []:
            for mot_coordinate in motif_positions:
                context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2])
                context.rectangle(mot_coordinate[0] + 101,
                                gene_ypos - 15,
                                motif_len,
                                30)
                context.fill()
        else:
            pass
        return context
    
    def draw_legend(self, context, motif_color: tuple, motif_name: str, motif_legend_ypos: int, canvas_width: int):
        '''This function draws the Legend on the canvas. Input the RGB value for the color of the motif as a tuple,
        the motif name, the y position for the motif entry in the legend and the width of the canvas.'''
        context.set_source_rgb(0,0,0)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(14)
        context.move_to(canvas_width - 150, 30)
        context.show_text("Legend")
        context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2])
        context.rectangle(canvas_width - 150, motif_legend_ypos, 15, 15)
        context.fill()
        context.set_source_rgb(0,0,0)
        context.move_to(canvas_width - 130, motif_legend_ypos + 10)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(14)
        context.show_text(motif_name)
        return context

# define functions
def oneline_fasta(multi_line_fa: str) -> dict:
    '''This function takes a fasta text file as an argument and creates a new fasta file with duplicate
    information but multi-line sequences are converted to one line such that format alternates between 
    header line and sequence line.'''
    seq = ''
    header = ''
    seq_dict = {}
    with open(multi_line_fa, 'r') as gh:
        for line in gh:
            line = line.strip('\n')
            if line.startswith('>'):
                if seq == '':
                    header = line
                else:
                    seq_dict[header] = seq
                    seq = ''
                    header = line
            else:
                seq += line
        seq_dict[header] = seq
    return seq_dict

def motifs(motifs_txt_file: str) -> list:
    '''This function reads through a text file containing a list of motifs and returns the motifs as a list of 
    strings.'''
    motifs_list = []
    with open(motifs_txt_file, 'r') as mf:
        for line in mf:
            line = line.strip('\n')
            motifs_list.append(line)
    return motifs_list

def calc_canv_width(fasta_dictionary: dict) -> int:
    '''Calculates the width of the canvas according to the base pair length of the longest gene sequence.'''
    for key in fasta_dictionary:
        gene_leng = []
        gene_leng.append(len(fasta_dict[key]))
    max_len = max(gene_leng) + 300
    return max_len

def calc_canv_height(fasta_dictionary: dict) -> int:
    '''Calculates the height of the canvas according to the number of genes in the fasta file.'''
    gene_count = 0
    for key in fasta_dictionary:
        gene_count += 1
    max_genes = gene_count*100 + 200
    return max_genes

def get_png_name(fasta_file: str) -> str:
    '''Gets the name of the fasta file, which will become the name of the png of the diagram.'''
    output_png = re.findall("Figure.+.fasta", fasta_file)
    return output_png[0][:-6]

# the algorithm
# read motifs text file in as a list
mot_list = motifs(motifs_file)
# read fasta file in as a dictionary (key: header, value: seq)
fasta_dict = oneline_fasta(fasta_file)
# parse output file name from fasta file name 
output_png = get_png_name(fasta_file)

# initialize y positions for gene name and diagram so they increment appropriately 
name_ypos, diagram_ypos = -20, 0
# initialize list of colors that will be assigned to each motif
colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (170, 0, 255), (0, 200, 250)]
# calculate dimensions of canvas 
canv_width = calc_canv_width(fasta_dict)
canv_height = calc_canv_height(fasta_dict)

# set up canvas object
canvas = Draw(100, 0, canv_width, canv_height)

# create the surface
surface, context = canvas.init_canvas(canvas.init_xpos, canvas.init_ypos, canvas.canv_wid, canvas.canv_hgt)

# iterate through each gene in the fasta dictionary
for key in fasta_dict:
    name_ypos += 100
    diagram_ypos += 100
    seq = fasta_dict[key]
    gene_len = len(seq)
    # isolte the gene name from the header 
    gene_name = re.findall(">[A-Z0-9]+", key)
    gene_name = gene_name[0][1:]

    # create gene object
    gene1 = Gene(gene_name, gene_len, seq, diagram_ypos, name_ypos)

    # write gene name on canvas
    context = canvas.write_gene_name(context, gene1.name, canvas.init_xpos, gene1.name_ypos)
    # draw gene on canvas
    context = canvas.draw_gene(context, canvas.init_xpos, gene1.diag_ypos, gene1.len)
    # get positions of exons 
    exon_pos = gene1.exon_finder(seq)
    # draw black rectangles as exon positions
    context = canvas.draw_exon(context, exon_pos, gene1.diag_ypos)

    # initialize legend y position in preparation for incrementation in next for loop
    legend_mot_ypos = 20
    for i,item in enumerate(mot_list):
        mot_seq = item
        mot_col = colors[i]
        mot_len = len(item)
        legend_mot_ypos += 20

        # create motif object 
        motif1 = Motif(mot_seq, mot_col, mot_len, legend_mot_ypos)

        # get all possible sequences for the motif containing ambiguous bases 
        motif1_all = motif1.all_possible_seq(mot_seq)
        # get positions of all instances of the motif in the gene sequence
        motif_pos = motif1.unlisting_fxn(motif1.all_motif_positions(motif1_all, gene1.seq))
        # draw the motif 
        context = canvas.draw_motif(context, motif_pos, motif1.color, gene1.diag_ypos, motif1.length)
        # draw the corresponding entry on the legend for this current motif
        context = canvas.draw_legend(context, motif1.color, motif1.seq, motif1.legend_ypos, canvas.canv_wid)

# convert the canvas into a png image file
surface.write_to_png(output_png + ".png")