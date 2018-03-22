
© 2018 Merly Escalona (<merlyescalona@uvigo.es>)

University of Vigo, Spain, http://darwin.uvigo.es

[![Build Status](https://travis-ci.org/merlyescalona/simphycompress.svg?branch=master)](https://travis-ci.org/merlyescalona/simphycompress)

# SimPhy compress dataset

This program allows to compress the number of files and file sizes of a SimPhy run.

# Assumptions

- We are working under a [SimPhy](https://github.com/adamallo/simphy) Meaning, it follows hierarchical  [SimPhy](https://github.com/adamallo/simphy) 's folder structure and sequence
labeling.

To know more about the simulation pipeline scenario go to:

- [SimPhy: A comprehensive simulator of gene family evolution ](https://github.com/adamallo/simphy)

# Input

- [SimPhy](https://github.com/adamallo/simphy) folder path
- prefix of the existing [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files
- (optional) length of the N sequence that will be used to separate the sequences when concatenated

# Output

- Modifications are made INPLACE. Meaning, files are concatenated and gzipped in the same SimPhy folder. And so,
the other files are removed.

# Install

- Clone this repository

```
git clone git@github.com:merlyescalona/simphycompress.git
```

- Chance your current directory to the downloaded folder:

```
cd simphycompress
```

- Install:

```
python setup.py install --user
```

# Usage

Required arguments:
- `-s <path>, --simphy-path <path>`: Path of the SimPhy folder.
- `-ip <input_prefix>, --input-prefix <input_prefix>`: Prefix of the FASTA filenames.

Optional arguments:
- `-n <N_seq_size>, --nsize <N_seq_size>`
    Number of N's that will be introduced to separate the sequences selected. If the parameter is not set, the output file per replicate will be a multiple alignment sequence file, otherwise, the output will be a single sequence file per replicate consisting of a concatenation of the reference sequences selected separated with as many N's as set for this parameter.
- `-l <log_level>, --log <log_level>`
    Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:['DEBUG', 'INFO', 'WARNING', 'ERROR']. Default: 'INFO'.

Information arguments:
  - `-v, --version`: Show program's version number and exit
  - `-h, --help`:    Show this help message and exit
