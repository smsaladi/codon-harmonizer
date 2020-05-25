[![Build Status](https://travis-ci.com/smsaladi/codon-harmonizer.svg?branch=master)](https://travis-ci.com/smsaladi/codon-slim)
[![PyPI version](https://badge.fury.io/py/codon-harmonizer.svg)](https://badge.fury.io/py/codon-harmonizer)
![PyPI - Downloads](https://img.shields.io/pypi/dm/codon-harmonizer)
[![DOI](https://data.caltech.edu/badge/209835343.svg)](https://data.caltech.edu/badge/latestdoi/266429793)

codon harmonizer
================

A program that attempts to recode a gene by considering the relative usage
of each codon in it's host and selecting a codon with the nearest relative
usage in a target organism.

In typical codon optimizers, each codon of a gene of interest is converted to
the "best" codon for a target organism. Yet, wild-type sequences don't always
use the "best" codon in their host organism. This code adjusts for this by
selecting the codon of a target organism that most closely approximates the 
codon's usage in a source organism.

This code was intended as a fork of Bart Nijsse's
[`codon harmonizer`](https://gitlab.com/wurssb/Codonharmonizer)
but ended up a whole-scale rewrite.

Where referenced in academic work, you may cite this repository and may also
consider referencing [manuscript](doi.org/10.1371/journal.pone.0184355)
discussing Nijsse's work.

### Installation
```bash
pip install git+git://github.com/smsaladi/codonharmonizer.git
```

### Usage

Count codon usage in source and target organisms
```bash
wget -o Gvio.cds.fna https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/385/GCF_000011385.1_ASM1138v1/GCF_000011385.1_ASM1138v1_cds_from_genomic.fna.gz
wget -o Ecol_MG1655.cds.fna https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz

codonharmonizer Gvio.cds.fna --write_freqs > Gvio.freq.csv
codonharmonizer Ecol_MG1655.cds.fna --write_freqs > Ecol_MG1655.freq.csv
```

Use these reference sets to recode genes of interest
```bash
codonharmonizer test/example_gene.fasta --target Ecol_MG1655.freq.csv --source Gvio.freq.csv
```
