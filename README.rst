===================
MLST
===================

This project documents MLST service


Documentation
=============

## What is it?

The MLST service contains one perl script *MLST-1.8.pl* which is the script of the lates
version of the MLST service. The method enables investigators to determine the ST based on WGS data.

## Usage

To use the service servral data needs to be pre-installed.
*database*
*blast-2.2.26*
*Makefile*

The folder *database* includes all the MLST schemes and needs to be updataed to get the best results.
The datasets are extracted from the http://pubmlst.org/ webside weekly and can be downloaded from
http://cge.cbs.dtu.dk/services/data.php. The folder also includes the *mlst_schemes* file which is
used by the program but also contains the different MLST schemes provided and what to write
in the commandline to use a particular shceme. 

The folder *blast-2.2.26* includes blastall and formatdb which used by the *MLST-1.8.pl* script

The file *Makefile* installs the nesserary perl modules to run the *MLST-1.8.pl* script.

The program can be invoked with the -h option to get help and more information of the service.

```bash
Usage: perl MLST-1.8.pl [options]

Options:

    -h HELP
                    Prints a message with options and information to the screen
    -d DATABSE
                    The path to where you have located the database folder
    -blast BLAST
                    The path to the location of blast-2.2.26
    -i INFILE
                    Your input file which needs to be preassembled partial or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places. If not specified the program will
                    create a folder named 'output' in which the result files will be stored.
    -s SPECIES
                    The MLST scheme you want to use. The options can be found in the *mlst_schemes* file
                    in the *database* folder
```

## Example of use with the *database* and *blast-2.2.26* folder loacted in the current directory 
    
    perl MLST-1.8.pl -i INFILE.fasta -o OUTFOLDER -s ecoli

## Example of use with the *database* and *blast-2.2.26* folder loacted in antoher directory

    perl MLST-1.8.pl -d path/to/database -blast path/to/blast-2.2.26 -i INFILE.fasta -o OUTFOLDER -s ecoli
    

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/MLST/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/mlst/overview

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/mlst/overview.

Installation
=======

The scripts are self contained. You just have to copy them to where they should
be used. Only the *database* folder needs to be updataed mannually. 

Remember to add the program to your system path if you want to be able to invoke the program without calling the full path.
If you don't do that you have to write the full path to the program when using it.

Citation
=======

When using the method please cite:

Multilocus Sequence Typing of Total Genome Sequenced Bacteria.
Larsen MV, Cosentino S, Rasmussen S, Friis C, Hasman H, Marvig RL,
Jelsbak L, Sicheritz-Pont�n T, Ussery DW, Aarestrup FM and Lund O.
J. Clin. Micobiol. 2012. 50(4): 1355-1361.
PMID: 22238442         doi: 10.1128/JCM.06094-11
[Epub ahead of print]


License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
