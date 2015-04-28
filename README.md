#### Development: ![Development](https://magnum.travis-ci.com/LeeBergstrand/BioMagick.svg?token=gFspxRhLX7xmypaS1yi5&branch=develop) 
#### Master: ![master](https://magnum.travis-ci.com/LeeBergstrand/BioMagick.svg?token=gFspxRhLX7xmypaS1yi5&branch=master)

---------------------------------------------------------------------
BioMagick
=========
### A next generation bioinformatics file format converter.

BioMagick is a next generation bioinformatics file format converter built from Biopython's SeqIO, AlignIO and Phylo classes. It currently contain a class for auto-recognition of text and binary Bioinformatics file formats. The command line user interface (CLI) allows piped input from other UNIX programs and the input of multiple files. Automated tests have been created to validate that BioMagick is providing the correct output. Documentation will be in the form of a Github wiki.

Features
--------
- Capable of automatic file type recognition of one or more bioinformatic file formats or the bioinformatic file type of a text stream from stdin via BioID.

- Once input file type is identified Biomagick converts file or text stream to the user indicated file type (compatibility is determined by Biopython) and writes this converted file to disk or stdout (stdout can be used only if there was one input file).

- Capable of converting multiple different input files of different types to a single file type (Many to one conversion).

- Has multiprocessing support. By default creates one subprocess per CPU core and converts on input file per core in a work-crew manor.

Supported Operating Systems
---------------------------
Biomagick is actively devleloped and tested on:

- Mac OSX
- Linux (Primarily Ubuntu)
- Windows 7

Other unix operating systems and linux distros should work assuming that one can install Python and the required dependencies. 


Dependencies
------------
- **Python** 2.6, 2.7, 3.3, 3.4 or newer
- Python Packages:
	- **biopython** 1.63
	- **binaryornot** 0.3.0
	- **PyYAML** 3.11

### To install dependancies on UNIX operating systems:

1. Use a package manager to install python if it is not present in your OS.

	```
	sudo apt-get install python # Ex. for Debian based linux 
	```
2. Intall pip

	```
	python get-pip.py
	```
3. Change directory to the BioMagick directory and intall python dependencies via pip

	```
	cd ~/the_file_path/Biomagick/
	pip install -r requirements.txt
	```
	
### To install dependancies on Windows:
- TBA

Command Line Interface Overview
-------------------------------
    usage: BioMagick.py [-h] [-i INPATH [INPATH ...]] [-s] [-o OUTPATH]
                    [-f OUTFORMAT [OUTFORMAT ...]] [-a ALPHA] [-j JOBS]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPATH [INPATH ...], --input INPATH [INPATH ...]
                        A list of input file paths. If not specified, input is
                        read from stdin.
      -s, --stdout          Output result of single-file conversion to stdout.
      -o OUTPATH, --outdir OUTPATH
                        An output directory for output files. If not
                        specified, output is piped to stdout.
      -f OUTFORMAT [OUTFORMAT ...], --outfmt OUTFORMAT [OUTFORMAT ...]
                        A List of output file formats.
      -a ALPHA, --alphabet ALPHA
                        The alphabet to use for conversion (ambigdna,
                        unambigdna, exdna, ambigrna, unambigrna, prot,
                        exprot).
      -j JOBS, --jobs JOBS  The number of processes to use for multiple files
                        (defaults to the number of processor cores).
                        
Documentation
-------------
Our documentation can be found in the **[Wiki](http://github.com/LeeBergstrand/BioMagick/wiki)**.


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

BioMagick is open-sourced and released under [MIT License](http://en.wikipedia.org/wiki/MIT_License).

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