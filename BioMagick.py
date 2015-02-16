#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A next generation bioinformatics file format converter and sequence feature extractor.
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
# ============================================================================================================

from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
# from BioID import BioID
import click


@click.command()
@click.option('-in', '--input_files', help='A comma separated list of input file paths. If not specified string is read from stdin and decode using BioID')
@click.option('-outfmt', '--output_formats', help='A comma separated list of output file formats. If more than one output format is specified, one output file per output format will be created per input file.')
@click.option('-outdir', '--output_directory', type=click.Path(), help='Output directory for output files. If not specified pipes output to stdout.')
@click.option('-lang', '--language', default='DNA', help='Output file language as DNA, RNA or Protein.')
@click.option('-r', '--reverse', help='Output file will contain reverse-complement of input sequence. Compatible with only DNA and RNA input languages.')
@click.option('-feat', 'features', help='Comma separated list of features to be extracted from input format (Genbank Only)')
@click.option('-afeat', '--anti_features', help='Comma separated list of features to not be extracted from input format (Genbank Only)')
@click.option('-mark', '--makers', help='Comma separated list of REGEX used to extract features with a given feature name (protein name etc). (Genbank Only)')
@click.option('-amark', '--anti_makers', help='Comma separated list of REGEX used to not extract features with a given feature name (protein name etc). (Genbank Only)')
@click.option('-seqrange', help='Output sequences will be trimmed to serquence range. If extracting features, each extracted feature will be trimmed to sequence range')
@click.option('-featrange', help='Extract features from parent sequence over a given range of parent sequence.  (Genbank Only)')
@click.option('-g', '--group', help='Group data from input files into single output file. (aggregate)')
@click.option('-degap', '--group', help='Remove gap symbols from sequence alignment files.')
def cli(input_files):
	"""BioMagick is a next generation bioinformatics file format converter and sequence feature extractor built from Biopython's SeqIO, AlignIO and Phylo classes."""
	click.echo(input_files)

if __name__ == '__main__':
	cli()