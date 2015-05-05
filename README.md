#### Development: ![Development](https://magnum.travis-ci.com/LeeBergstrand/BioMagick.svg?token=gFspxRhLX7xmypaS1yi5&branch=develop) 
#### Master: ![master](https://magnum.travis-ci.com/LeeBergstrand/BioMagick.svg?token=gFspxRhLX7xmypaS1yi5&branch=master)

---------------------------------------------------------------------
BioMagick
=========
### A next generation bioinformatics file format converter.

BioMagick is a next generation bioinformatics file format converter built from Biopython's SeqIO, AlignIO and Phylo classes. It currently contains a class for auto-recognition of text and binary Bioinformatics file formats. The command line user interface (CLI) allows piped input from other UNIX programs and the input of multiple files. Automated tests have been created to validate that BioMagick is providing a correct output. Documentation can be found in the Github wiki.

Features
--------
- Capable of identifying the file type of both text and binary bioinformatic files as well as the file type of a text stream of a bioinformatic file piped into BioMagick via standard in.

- Capable of converting multiple input files of different types or a singular text stream into one or more output formats (one output file per output format per input file or text stream) using a single command.

- Capable of outputting the output of a single format conversion of a single input file or text stream to standard out.

- Capable of multiprocessing with up to one input file being converted per CPU core.

### The following formats are supported by BioMagick:

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

Supported Operating Systems
---------------------------
Biomagick is actively developed and tested on:

- Mac OSX
- Linux (Tested On Ubuntu 12.04 LTS Server Edition 64 bit)
- Windows 7

Other unix operating systems and linux distributions should work assuming that one can install Python and the required dependencies. 

Dependencies
------------
- **Python** 2.6, 2.7, 3.3, 3.4 or newer
- Python Packages:
	- **biopython** 1.63
	- **binaryornot** 0.3.0
	- **PyYAML** 3.11

### To install dependancies on UNIX operating systems:

1. Use a package manager to install python if it is not included in your OS.

	```
	sudo apt-get install python # Ex. for debian-based linux distributions
	```
2. Install pip

	```
	python get-pip.py
	```
3. Change directory to the BioMagick directory and install python dependencies via pip

	```
	cd ~/the_file_path/Biomagick/
	pip install -r requirements.txt
	```
	
### To install dependancies on Windows:
- TBA

Command Line Interface Overview
-------------------------------

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
Documentation
-------------
Our documentation can be found in the **[Wiki](http://github.com/LeeBergstrand/BioMagick/wiki)**. There is also an included user guide.

Trouble Shooting
----------------
If you forget the command line arguments try typing `python biomagick.py -h`

Please report issues to the **[issues page](http://github.com/LeeBergstrand/BioMagick/issues)**.  


File Manifest
-------------
	.
	├── BioID.py *
	├── BioMagick.py
	├── BioIDFormatInfo.yml *
	├── BioMagickFormatInfo.yml
	├── LICENSE.txt
	├── README.md
	├── __init__.py
	├── media/
	│   ├── New_Data_Flow.png
	│   ├── Old_Data_Flow.png
	│   ├── Parallelism.png
	│   ├── UML_Complete.png
	│   ├── UML_Current.png
	│   └── UML_Future.png
	├── reference/
	│   ├── Project_Proposal.md
	│   ├── Statistical_Learning_for_File_Type_Identification.pdf
	│   ├── fileExtentions/
	│   │   ├── AlignIO.csv
	│   │   ├── Phylo.csv
	│   │   └── SeqIO.csv
	│   └── unix_prog_design.pdf
	├── requirements.txt
	└── testing/
	    ├── conversion_tests.yml
	    ├── format_tests.csv *
	    ├── testFiles/ **
	    ├── test_BioID.py *
	    └── test_BioMagick.py
	    
	* These files may be moved to their own repository if BioID is released as its own package.
	** Contains many test files

Licence
-------

BioMagick is open-source and released under [MIT License](http://en.wikipedia.org/wiki/MIT_License).

	The MIT License (MIT)
	
	Copyright (c) 2014 Lee Bergstrand and Matt McInnes
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.

Authors
-------
[Lee Bergstrand](http://github.com/LeeBergstrand)

[Matt McInnes](https://github.com/Krailon)
