## PROGRAM USED FOR THE STUDY OF TROPHIC MODES IN PROTISTS

## Prerequisites

The use of [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) is strongly recommended.

### Miniconda Installation

> *Windows:* https://conda.io/projects/conda/en/latest/user-guide/install/windows.html

> *MacOS:* https://conda.io/projects/conda/en/latest/user-guide/install/macos.html

> *Linux:* https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Diamond Installation

This program uses [Diamond](https://github.com/bbuchfink/diamond), you can install it *via* conda:

```
conda install -c bioconda diamond
```

## Download the program

1. Clone git repository

> HTTPS link

```
git clone https://github.com/Dylkln/EGmnhn.git
```

> SSH link

```
git clone git@github.com:Dylkln/EGmnhn.git
```

2. Initialize conda environment using the file *environment.yml*

```
conda env create --file environment.yml
```

3. Activate conda environment

```
conda activate trophomics
```

## Usage

This is a workflow using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

### Configuration

You must complete the config file ***config.yml*** to launch every step you need.
Leave a blank for any steps you don't want.


### Utilisation

You can familiarize yourself with the program using the test script.
This script is located on EGmnhn/test. You must have at least 6GB of RAM to run the test.

To use it, you can launch this command :

```
snakemake -c <NUMBER OF CORES>
```

OR

```
snakemake -c all
```

The first command asks snakemake to use a specific number of cores, the second asks to use all cores
available

If you use this command in the test repository, it will treat all test files (transdecoder and
annotation), will filter the network keeping sequences with at least 80% overlap and with at least
80% identity. The column names used as attributes are "peptides", "length", and "database".

### Data used by the program

All you data must be standardized, please refer to the config.yml.

Your Input must be a path to a fasta file or a file containing the list of path to fasta files.
Annotation files must contain in the first column all **ID/Name** of sequences contained in fasta
format files.
The test data are part of a bigger dataset called [MetDB](http://metdb.sb-roscoff.fr/metdb/).