
© 2018 Merly Escalona (<merlyescalona@uvigo.es>)

University of Vigo, Spain, http://darwin.uvigo.es


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
