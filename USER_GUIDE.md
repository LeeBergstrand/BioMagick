# BioMagick User Guide

## Comand Line Interface Overview

	usage: BioMagick.py [-h] [-i INPATH] [-s] [-o OUTPATH] [-f OUTFORMAT]
	                    [-a ALPHA] [-j JOBS]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  
	  -i INPATH, --input INPATH
	                        A comma-separated list of input file paths. If not
	                        specified, input is read from stdin.
	  
	  -s, --stdout          Output result of single-file conversion to stdout.
	  
	  -o OUTPATH, --outdir OUTPATH
	                        An output directory for output files. If not
	                        specified, the current working directory is used.
	  
	  -f OUTFORMAT, --outfmt OUTFORMAT
	                        A comma-separated list of output file formats.
	  
	  -a ALPHA, --alphabet ALPHA
	                        The alphabet to use for conversion (ambigdna,
	                        unambigdna, exdna, ambigrna, unambigrna, prot,
	                        exprot).
	  
	  -j JOBS, --jobs JOBS  The number of processes to use for multiple files
	                        (defaults to the number of processor cores).
                        
## Basic Usage


BioMagick requires, at the very least, two things to function:

* input file(s)
* output format(s)
 
### Input

Input can be specified in one of two ways. Single input files can be specified via the **-i/--input** argument as seen below.

```
./BioMagick.py -i path/to/file1 ...
```
 
Alternatively, in the case of multiple input files, one can specify each input file path in a comma-separated list.

```
./BioMagick.py -i path/to/file1,path/to/file2,path/to/file3 ...
```

For specifying single files another option is also available. The desired input file's contents can be piped directly into BioMagick for conversion as shown below.

```
sed 's/replacee/replacer' input_file | ./BioMagick.py ...
```

### Output


Output is controlled by three arguments in BioMagick. The first is **-o/--outdir**, which specifies what directory to save the output file(s) in. This argument is optional, and will default to the current working directory if omitted.

```
./BioMagick.py -i infile1,infile2,infile3 -o /home/someuser/datafiles ...
```

The second output argument is **-f/--outfmt**, which specifies the output format(s) to convert to in a comma-separated list. The output format name must be specified in lowercase, as shown below.

```
./BioMagick.py -i infile1,infile2,infile3 -f fasta,stockholm ...
```
**The following formats are supported by BioMagick:**

| Format Name | BioMagick Name | Description                                                                                                                    |
|-------------|----------------|--------------------------------------------------------------------------------------------------------------------------------|
| FASTA       | fasta          | This refers to the input FASTA file format introduced for Bill Pearson's FASTA tool, where each record starts with a ">" line. |
| FASTQ       | fastq          | FASTQ files are a bit like FASTA files but also include sequencing qualities.                                                  |
| PHYLIP      | phylip         | An alignment format used by the Phylip software package.                                                                                                           |
| Clustal X   | clustal        | The alignment format of the Clustal X aligner.                                                                                  |
| Genbank     | genbank        | The GenBank or GenPept flat file format.                                                                                       |
| EMBL        | embl           | The EMBL flat file format.                                                                                                     |
| SeqXML      | seqxml         | An XML-based sequence file format.                                                                                              |
| PhyloXML    | phyloxml       | An XML-based phylogenetic tree file format.                                                                                      |
| NeXML       | nexml          | An XML version of the Nexus format.                                                                                            |
| Nexus       | nexus          | The NEXUS multiple alignment format, also known as PAUP format.                                                                |
| Newick      | newick         | An phylogenetic tree file format used by the Newick software package.                                                           |
| Stockholm   | stockholm      | The Stockholm alignment format is also known as PFAM format.                                                                   |
| SFF         | sff            | Standard Flowgram Format (SFF) binary files produced by Roche 454 and IonTorrent/IonProton sequencing machines.                |

The third output argument is **-s/--stdout**, which specifies that the output should be printed to the standard output stream after conversion. This option only works with single-file conversions since there's no easy way to differentiate between multiple output files in a stream.

```
./BioMagick.py -i infile -s | sed 's/replacee/replacer' ...
```

Using piped input and output, BioMagick can be used for complex stream processing. For example the below shell command uses bedtools to convert a BAM file to FASTQ. The bedtools output is piped into BioMagick, converted to FASTA, altered via SED and written to a file.    

```
bedtools bamtofastq  -i sequences.bam -fq - | ./BioMagick.py -f fasta ... | sed 's/replacee/replacer' > sequences.fasta
```

### SeqIO Alphabet


The alphabet argument **-a/--alphabet** is used to specify the input sequence alphabet. This is only needed when converting to a format that uses BioPython's **SeqIO** conversion class. Any conversions using **AlignIO** or **Phylo** will ignore any specified alphabet.

```
./BioMagick.py -i infile1,infile2,infile3 -f fastq -a unambigdna ...
```

### Multiprocess Job Control


The job argument **-j/--job** allows control over how many parallel conversion processes are used when converting multiple files. By default, the number of CPU cores is used if none is specified.

```
./BioMagick.py -i infile1,infile2,infile3,infile4,infile5,infile6 -f nexus -j 4 ...
```
