To use the tool online go to:<br />
http://codonharmonizer.systemsbiology.nl

Please cite this tool as:<br />
Claassens NJ, Siliakus MF, Nijsse B, Spaans SK, Creutzburg SCA, Schaap PJ, et al. Improving heterologous membrane protein production in Escherichia coli by combining transcriptional tuning and codon usage algorithms. PLoS One. 2017

Requirements: Python 3

```
Harmonize genes for a target organism.

usage:
codonharm.py -f <(multi)fasta_file> -o <output file> -s <frequency_file> -t <frequency_file>,<frequency_file>, etc..

Harmonize your genes for a target organism. See codonfrequencies_from_cds.py to generate frequency files.

optional arguments:
  -h, --help            show this help message and exit

Input:
  -f FASTA, --fasta FASTA
                        DNA (multi)fasta file
  -s SOURCE, --source SOURCE
                        Source Organism eg. Eco_MG1655
  -t TARGETS, --target TARGETS
                        Target organism(s) eg. Eco_MG1655. Can be a comma separated list.
  -o NAME, --output NAME
                        Output filename (.zip)


usage:
codonfrequencies_from_cds.py -n <name> -o <filename> <CDS-fasta>

Generate a frequency file from a CDS fasta file used for the codonharmonizer

positional arguments:
  CDS-FASTA          DNA multi-fasta file of protein coding genes

optional arguments:
  -h, --help         show this help message and exit
  -n, --name NAME    Name of the organism
  -o, --output FILE  output file (.csv)
  -q, --quiet        Ignore warnings
  
  
  
```