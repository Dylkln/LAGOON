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
conda activate trophic
```

## Usage

### Arguments

```
python egmnhn.py -tr <TRANSDECODER_FILES> -an <ANNOTATION_FILES> -ov <OVERLAP> -id <IDENTITY> -cn <COLUMN_NAMES>
```

All arguments are optional.

- ***TRANSDECODER_FILES*** : either a path to find files in fasta format or a file containing the list of all path
to files in fasta format.
- ***ANNOTATION_FILES*** : either a path to find annotation files or a file containing the list of all path
to annotation files.
- ***OVERLAP*** : The percentage of overlap between sequences used as filtration parameter.
- ***IDENTITY*** : The percentage of identity between sequences used as filtration parameter.
- ***COLUMN_NAMES*** : The names of the columns wanted as attributes for the sequence similarity network.

### Example of utilisation

You can familiarize yourself with the program using the test script.
This script is located on EGmnhn/test.

To use it, you can launch this command :

```
python test.py -tr ./tr_files_test/* -an ./an_files_test/* -ov 80 -id 80 -cn peptides length database
```

This command will treat all test files (transdecoder and annotation), will filter the network keeping sequences with at
least 80% overlap and with at least 80% identity. The column names used as attributes are "peptides", "length", and
"database".

### Data used by the program

All you data must be standardized. All files used in *<TRANSDECODER_FILES>* argument must be in fasta format.
Annotation files must contain in the first column all **ID/Name** of sequences contained in fasta format files.
The test data are part of a bigger dataset called [MetDB](http://metdb.sb-roscoff.fr/metdb/).