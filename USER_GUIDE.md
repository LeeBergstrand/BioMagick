# BioMagick User Guide


```

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

```


## Basic Usage


BioMagick requires, at the very least, two things to function:

* input file(s)

* output format(s)


### Input


Input can be specified in one of two ways. For multiple files, the **-i/--input** argument

is used to specify each input file path in a comma-separated list as shown below.


```shell

./BioMagick.py -i path/to/file1,path/to/file2,path/to/file3 ...

```


For specifying single files, the **-i**/**--input** argument can still be used, but another option 

is also available that enables I/O chaining. The desired input file's contents can be 

piped directly into BioMagick for conversion as shown below.


```shell

sed 's/replacee/replacer' input_file | ./BioMagick.py ...

```


### Output


Output is controlled by three arguments in BioMagick. The first is **-o/--outdir**, which specifies 

what directory to save the output file(s) in. This argument is optional, and will default to the current 

working directory if omitted.


```shell

./BioMagick.py -i infile1,infile2,infile3 -o /home/someuser/datafiles ...

```


The second output argument is **-f/--outfmt**, which specifies the output format(s) to convert to in a 

comma-separated list. The output format name must be specified in lowercase, as shown below.


```shell

./BioMagick.py -i infile1,infile2,infile3 -f fasta ...

```


The third output argument is **-s/--stdout**, which specifies that the output should be printed to

the stardand output stream after conversion. This option only works with single-file conversions since

there's no easy way to differentiate multiple output files in a stream and save them to different filenames.


This option is the logical opposite of the piped input and functions to allow insertion of BioMagick into

Bioinformatic pipelines as shown below.


```shell

./BioMagick.py -i infile -s | sed 's/replacee/replacer' ...

```


It can also be used to do successive conversions on the same file(s) by piping BioMagick's output into 

further instances of BioMagick as new input repeatedly as shown below.


```shell

./BioMagick.py -i infile1,infile2,infile3 -f genbank -s | ./BioMagick.py -f fastq -s | ./BioMagick.py -f fasta ...

```


### SeqIO Alphabet


The alphabet argument **-a/--alphabet** is used to specify a conversion alphabet. This is only needed when 

converting to a format that uses BioPython's **SeqIO** conversion class. Any conversions using **AlignIO** 

or **Phylo** will ignore any specified alphabet.


```shell

./BioMagick.py -i infile1,infile2,infile3 -f fastq -a unambigdna ...

```


### Multiprocess Job Control


The job argument **-j/--job** allows control over how many parallel conversion procesesses are used when

converting multiple files. By default, the number of CPU cores is used if none is specified.


```shell

./BioMagick.py -i infile1,infile2,infile3,infile4,infile5,infile6 -f nexus -j 4 ...

```
