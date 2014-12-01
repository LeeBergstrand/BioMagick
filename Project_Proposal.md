BioMagick
=========
### A next generation bioinformatics file format converter and sequence feature extractor.    


## Introduction

### Next Generation Sequencing

![Cost of sequencing in the past decade.](http://upload.wikimedia.org/wikipedia/commons/thumb/b/b7/DNA_Sequencing_Cost_per_Genome_Over_Time.jpg/380px-DNA_Sequencing_Cost_per_Genome_Over_Time.jpg)

With the democratization of [next generation sequencing (NGS)](http://en.wikipedia.org/wiki/DNA_sequencing#Next-generation_methods) technologies over the past decade the cost of sequencing an organism's genome is dropping drastically and therefore sequencing is an ever increasingly attractive analytical method for biologist everywhere. The amount of sequencing data being generated is also increasing due to NGS technologies outputting tens of gigabytes of data per run. To counter act this new torrent of data computational biologists and bioinformaticians have begun to develop bioinformatic pipelines.

### Bioinformatics Pipelines

![Example bioinformatics pipeline.](http://www.broadinstitute.org/~bhaas/euk_annot_pipeline.png)

Bioinformatic pipelines are designed to automate the processing, filtering and analyzing of sequencing data. They often consist of several command line tools pipelined together with shell scripts or build creation tools such as [GNU Make](http://en.wikipedia.org/wiki/Make_(software). In these pipelines the output from one bioinformatic tool is passed along to another, with each step in the pipeline acting to transform the sequence data into a more useful form. Often each tool in the pipeline was developed on its own and therefore creates inputs and outputs incompatible with other tools in the system. As a result "glue" code has to be developed in order to transform the output from one bioinformatic tool into an input compatible with another.

### Bioinformatics Libraries

Several open source bioinformatics libraries have been created over the past three decades in order to simplify the creation of custom bioinformatic software. The two most prominent examples are [Biopython](http://biopython.org/wiki/Main_Page) and the legacy [BioPerl](http://www.bioperl.org/wiki/Main_Page) although other libraries have been created or are being created such as [BioJava](http://biojava.org/wiki/Main_Page), [BioRuby](http://bioruby.org), and [Bionode](http://bionode.github.io/bionode-website/). There are also a variety of C/C++ libraries both new and old. Often classes from these libraries are incorporated into simple scripts that transform one bioinformatics format to another, extract certain sequence features (for example protein coding genes) or provide the computational biologist with simple stats about their sequence data. These small scripts act as the "glue code" in the larger bioinformatic pipeline.

#### BioMagick will be designed to replace these scripts with a universal sequence conversion program built on top of open source libraries and possessing an elegant [command line user interface (CLI)](http://en.wikipedia.org/wiki/Command-line_interface).                

## Motivation

##### Current Standard

The current standard for converting bioinformatics formats is writing your own scripts using classes from Biopython or BioPerl both of which can take care of file parsing and conversion. Unexperienced biologists often end up with one script per conversion and this results in them having to maintain hundreds of one-off scripts created by themselves or passed along to them by their co-workers. ([here is an example repository](https://github.com/carden24/Bioinformatics_scripts)). In addition, biologists often rush the creation of these scripts so they can move on to analyzing the biology. This often leads to conversion errors which may skew the researcher's results. 

As the amount of biological sequence data is exponentially increasing, computational biologists are starting to move towards automated bioinformatic pipelines. However, building these pipelines off of a series of custom conversion scripts (often with the exact same internal data flow) is ultimately unsustainable and as a result an integrated bioinformatic file format cross-conversion tool will be required in the future.  

##### Existing Programs

A variety of bioinformatic file format cross-conversion programs have been developed in the past however most have failed to obtain significant usage. These programs all have the following characteristics:

- A limited number of input/output file formats
- Can only take one input file at a time
- Unintuitive and/or non-[POSIX](http://en.wikipedia.org/wiki/POSIX) CLI
- Were built using "older" programming languages such as C, Java or shell in a unmaintainable way
- Do not utilize existing open-source libraries that are modern and maintained
- Generally do not auto-recognize the input file format  
 
All of the bioinformatic file format conversion programs that the author has identified have either been web applications which can only convert the contents of a single text file that has been copy-pasted into a web form and sent to a server for conversion (ex. [Sequence Format Converter](http://genome.nci.nih.gov/tools/reformat.html), [readseq-web](http://iubio.bio.indiana.edu/cgi-bin/readseq.cgi), [fmtseq](http://www.bioinformatics.org/JaMBW/1/2/index.html)) or command line programs which are no longer maintain or have lost availability (ex. [ReadSeq](http://sourceforge.net/projects/readseq/), [BIOFFORC](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639671/#!po=16.6667)). The most prominent and feature rich converter is [ReadSeq](http://sourceforge.net/projects/readseq/), which is a command line program. It is also utilized as the conversion program on the server side for a variety of bioinformatic file conversion web apps. ReadSeq is over a decade old as its last source code update was on June 4, 2004. ReadSeq was also written in a old version of Java that is no longer supported by modern operating systems. Overall, there are no conversion programs based off of modern Bioinformatic libraries and most have little to no documentation.

##### Inspiration

The primary inspiration for BioMagick is the [ImageMagick](http://www.imagemagick.org) command line bitmap image processing suite. In a similar manner to ImageMagick, BioMagick will facilitate the cross-conversion between a variety of bioinformatic file formats through a single command line interface.      

## Project Summary

BioMagick will be a next generation bioinformatics file format converter and sequence feature extractor built from Biopython's [SeqIO](http://biopython.org/wiki/SeqIO), [AlignIO](http://biopython.org/wiki/AlignIO) and [Phylo](http://biopython.org/wiki/Phylo) classes. It will contain a class for auto-recognition of text based Bioinformatics formats. It will also possess value added features such as the ability to extract sequence features such as protein coding or 16s rRNA genes and the ability to display certain statistics about an input file, for example how many DNA sequences are contained within. The command line user interface (CLI) will be POSIX/GNU standard and will allow piped input from other UNIX programs and the input of multiple files. Automated tests will be created to validate that BioMagick is providing the correct output. Documentation will be in the form of a Github wiki and a UNIX Man page. BioMagick will be open source and released under [MIT License](http://en.wikipedia.org/wiki/MIT_License).

## Project Details

### Architecture and Environment

##### Languages and Platforms

BioMagick will be designed to be run on [LINUX](http://en.wikipedia.org/wiki/Linux) and [Mac OSX](http://en.wikipedia.org/wiki/OS_X). The programming language we will use is [Python](https://www.python.org) and we will support the following Python interpreters:

- CPython 2.6, 2.7, 3.3, 3.4 -- see [http://www.python.org](http://www.python.org)

- PyPy 2.0, 2.1, 2.2, 2.3 -- see [http://www.pypy.org](http://www.pypy.org)

We will also support the Python 2.7 and Python 3.0-3.4 versions of the Python language. Do to Python being an interpreted language, BioMagick should run on most hardware platforms supported by Linux (x86, Power, ARM etc.). However, our testing will only guarantee full functionality on [x86](http://en.wikipedia.org/wiki/X86) hardware (Intel and AMD)

##### Generalized Architectural Diagram

![Generalized Architecture Diagram](https://github.com/LeeBergstrand/BioMagick/raw/master/media/DataFlow.png)

##### Biopython

A core component of BioMagick will be the [Biopython](http://biopython.org/wiki/Main_Page) library. Specifically, the [SeqIO](http://biopython.org/wiki/SeqIO), [AlignIO](http://biopython.org/wiki/AlignIO) and [Phylo](http://biopython.org/wiki/Phylo) classes from Biopython will be heavily utilized for parsing, converting and writing bioinformatic file formats. They can also be used for various other tasks  such as extracting protein coding genes.

Some useful Biopython features:

- A generic sequence record object for storing genomic sequences along with their features (for example the position of genes within the sequence).
- Parsers for parsing various bioinformatic file formats into sequence record objects.
- Writers for serializing sequence record objects into various bioinformatic file formats.

######  Biopython code example #1:

```
from Bio import SeqIO
handle = open("example.fasta", "rU") # Opens a file handle.	
	
# Parses input file which contains multiple FASTA formatted sequences.
# Parses each FASTA formatted sequence into a generic sequence record object.
for record in SeqIO.parse(handle, "fasta") : 
	print record.id # Prints the id of each sequence.
handle.close()
```
	
######  Biopython code example #2:

```
from Bio import SeqIO
input_handle = open("cor6_6.gb", "rU") 
output_handle = open("cor6_6.fasta", "w")
	
# Parses each Genbank file into a generic sequence record object.
sequences = SeqIO.parse(input_handle, "genbank") # Opens a file handle.
	
# Write sequence record objects to file in FASTA format.
SeqIO.write(sequences, output_handle, "fasta")
 
output_handle.close()
input_handle.close()
```

Biopython can also be used to extract sequence features from sequences if the input file has these features annotated as part of their file type. For example, when the Genbank file format, which contains annotations for certain sequence features, is parsed into a sequence record object these annotations are stored within the object and can be accessed and extracted by the programmer. 

######  Biopython code example #3:

``` 
for feature in seqRecord.features:
	# If feature is a Coding DNA Strand
	if feature.type == "CDS": 
		# Gets the feature location
		CDSLocal = feature.location
		# Gets sequence quantifiers. 
		featQualifers = feature.qualifiers 
		
		# Gets sequence quantifiers.
		gene = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
		product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
		proteinID = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
		locus = str(featQualifers.get('locus_tag','no_locus_tag')).strip('\'[]')
```

The author has used Biopython to create software systems in the past. Here are some examples:

- [Genbank-Downloaders](https://github.com/LeeBergstrand/Genbank-Downloaders) - A series of small Biopython scripts for downloading sequence data off NCBI's Genbank.
- [BackBLAST_Reciprocal_BLAST](https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST) -  A reciprocal BLAST program for filtering down BLAST results to best bidirectional hits. 
- [HMMER-DB](https://github.com/LeeBergstrand/HMMER-DB) - A program for searching, extracting and storing proteins that match specific Hidden Markov Models created by HMMER.

##### Automatic Identification of Bioinformatic File Formats.

One of the features the author seeks to include into BioMagick is the ability to automatically recognize textual bioinformatic file formats. Once identified the file type of an input file would be passed on to Biopython. Automatic identification of file format was not included in Biopython as its author's made the conscious decision that input file format has to be explicitly stated before parsing. 

[BioPerl's SeqIO](http://bioperl.org/wiki/HOWTO:SeqIO) class (the Perl equivalent to Biopython's SeqIO class) is capable of automatically identifying a limited set of input file types. For identification of file type one could reverse engineer the BioPerl SeqIO class and identify the methods it uses for identifying file type or perhaps we could take a more generic approach using a library like [libmagic](http://linux.die.net/man/3/libmagic) to identify the input file.

##### Command Line Interface Libraries

A command line library will be used to parse user input. Some possible libraries are [Clint](https://pypi.python.org/pypi/clint/), [Click](http://click.pocoo.org/3/) or [Cliff](http://cliff.readthedocs.org/en/latest/).

##### Testing and Continous Integration

[Automated regression tests](http://en.wikipedia.org/wiki/Regression_testing) will be created to insure proper identification of bioinformatic file types. The [unittest](https://docs.python.org/2/library/unittest.html#module-unittest) module, which is built into Python's standard library, will be used as a testing framework. [Travis CI](https://travis-ci.org), a continuous integration service that is free for students and open source projects, will be used to run tests on each build of BioMagick. Biopython already includes automated regression testing for file parsing, writing and conversion by its SeqIO, AlignIO and Phylo classes, therefore we may not need to test this functionality. 

In later stages of the project we will have users beta test Biomagick. Potential beta testers are the [Mohn Lab](http://www.cmde.science.ubc.ca/mohn/) (the author's previous employer) at University of British Columbia and the [Van Hamme Lab](http://faculty.tru.ca/jvanhamme/) at Thompson Rivers University.

##### Version Control, Repository Hosting and Open Source Licensing

We will be using the [Git](http://git-scm.com) version control system for source control. Biomagick will be hosted in a private repository on [Github](https://github.com). In the later stages of testing we may move it to a public repository. BioMagick will be open-sourced and released under [MIT License](http://en.wikipedia.org/wiki/MIT_License) by the end of the winter semester.

### Implementation Issues and Challenges

- We would not want to write individual parsers for hundreds of input file formats, however creating a module for identifying multiple text file formats in a high throughput manor would be challenging. 

- Writing test cases will be time consuming. We should develop some abstractions to help us write these test cases in a more universal manor.

- Some Bioinformatics formats such as [Genbank](http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) contain significantly more sequence data (for example information about the genes contained within a parent sequence) than other, more raw formats such as [FASTA](http://en.wikipedia.org/wiki/FASTA_format) (which only contains raw sequences). It will be important to decide which sequence formats can be effectively converted into another. For example, Biopython is capable of "up-converting" low information file formats into high information file formats, for example FASTA to Genbank, however we should give the user a warning when this occurs to notify the user that the output file will not contain full information.

- There are two versions of the Python language available 2.x and 3.x. Python 3.0 was released in 2008 and featured some major changes that broke backwards compatibility with older versions of Python. Python 2.x is no longer being maintained. Therefore, Python 2.x is legacy and Python 3.x is the present and future of the language. Python 2.x is still heavily used in production systems and has a larger number of available external packages. There are ways to write universal libraries that are compatible with both families of Python and this is the path Biopython has followed. We also intend to write BioMagick in the same manor though this may have some unforeseen challenges. 

### Deliverables

- 4 + 1 Views Software Architecture Model
- Test plan and tests:
	- Unit tests for the bioinformatic file identification class.
	- Integration tests
	- Acceptance tests
- Table of bioinformatic file formats supported by Biopython and these files respective filename extensions.
- Prototype program using filename extensions for automatic file type identification.
- Fully developed bioinformatic file identification class.
- Release version of program
- Documentation:
	- For the file identification class.
	- For the command line interface. 
- Poster describing work completed.

### Putative Timeline

| Task                                                                                                       | Week | Month    |
|------------------------------------------------------------------------------------------------------------|------|----------|
| 4 + 1 Views Software Architecture Model                                                                    | 3    | January  |
| Table of bioinformatic file formats supported by Biopython and these file type's respective filename extensions. | 4    | January  |
| Test plan document                                                                                         | 1    | February |
| Unit tests for the bioinformatic file identification class                                                 | 2    | February |
| Integration and Acceptance tests                                                                           | 3    | February |
| Prototype program using filename extensions for automatic file type identification                         | 1    | March    |
| Fully developed bioinformatic file identification class                                                    | 3    | March    |
| Open sourcing the project and open beta testing                                                              | 4    | March    |
| Refactor and improve documentation                                                                         | 1    | April    |
| Review beta feedback and bug fix                                                                           | 2    | April    |
| Final Review Poster                                                                                        | 3    | April    |
| Final Release                                                                                              | 4    | April    |

### Student Skills and Experience Developed

- Experience building POSIX compliant command line programs
- Writing automated regression tests
- Running tests using Travis CI
- Utilizing technologies for identifying textual file formats
- Greater experience with Biopython
- Experience in creating an open source project
- Writing Python code that is compatible with both Python 2 and Python 3

## Conclusion

BioMagick will be a command line program for the conversion of various bioinformatics file formats. It will also have code for extracting certain sequence features from genomic data using Biopython. By building on existing expansive libraries such as Biopython we hope to reduce our time to completion thus allowing us to focus on other value added features such as automatic input file identification. BioMagick shall be completed in a directed studies course between January and April 2015.  
