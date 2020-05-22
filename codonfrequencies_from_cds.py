#!/usr/bin/python3

import sys
import argparse

gencode_11 = {
    "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
    "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
    "TTA": "L", "TCA": "S", "TAA": "-", "TGA": "-",
    "TTG": "L", "TCG": "S", "TAG": "-", "TGG": "W",
    "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
    "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
    "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
    "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"}


def parse_options():
    parser = argparse.ArgumentParser(
        description='Generate a frequency file from a CDS fasta file used for the codonharmonizer')

    parser.add_argument(dest="fasta_filepath",
                        help="DNA multi-fasta file of protein coding genes", metavar="CDS-FASTA")
    parser.add_argument("-n, --name", dest="name", required=True,
                        help="Name of the organism", metavar="NAME")
    parser.add_argument("-o, --output", dest="output",
                        required=True, help="output file (.csv)", metavar="FILE")
    parser.add_argument("-q, --quiet", dest="quiet",
                        action='store_true', help="Ignore warnings")

    inputs = parser.parse_args()

    return inputs


def get_sequences(fasta_contents):
    """Given a fasta file, return a dictionary of the entries"""
    try:
        sequence_dic = {}
        header = ""
        for line in fasta_contents:
            if line[0] == '>':
                header = line.strip()[1:]
                sequence_dic[header] = ""
            else:
                clean_seq = line.strip().upper().replace("U", "T")

                if header == "":
                    raise ValueError()

                sequence_dic[header] += clean_seq

        return sequence_dic

    except ValueError as err:
        sys.stderr.write("Not a valid DNA fasta file (missing header)")
        sys.exit()


def split_to_codons(sequence, header, quiet):
    """Returns a codon list from a sequence reading frame +1 """
    valid_bases = "ATCG"
    codon_sequence_list = []
    if len(sequence) % 3 == 0:
        codon_sequence_list = [sequence[i:i+3]
                               for i in range(0, len(sequence), 3)]
    elif not quiet:
        print("NOT USED: Partial sequence >"+header +
              " not divisible by complete codons")

    codon_sequence_list_clean = []
    for codon in codon_sequence_list:
        if all(char in valid_bases for char in codon):
            codon_sequence_list_clean.append(codon)
        elif not quiet:
            print("Removed codon sequence", codon,
                  "containing non DNA letter from:", header)
    return codon_sequence_list_clean


def main(argv):
    inputs = parse_options()

    try:
        fastafile = open(inputs.fasta_filepath, "r").readlines()
    except (OSError, IOError) as err:
        print("Unable to open input fasta file")
        sys.exit()

    fasta = get_sequences(fastafile)
    if len(fasta) < 100:
        print("WARNING: Number of genes < 100")
    frequentie_file = open(inputs.output, "w")
    frequentie_file.write("NAME,"+inputs.name)

    codon_frequencies = {}
    for entry in fasta:
        codons = split_to_codons(fasta[entry], entry, inputs.quiet)
        for codon in codons:
            if codon not in codon_frequencies.keys():
                codon_frequencies[codon] = 1
            else:
                codon_frequencies[codon] += 1

    for codon in codon_frequencies:
        AA = gencode_11[codon]
        frequentie_file.write("\n"+codon+","+AA+"," +
                              str(codon_frequencies[codon]))

    frequentie_file.close()


if __name__ == "__main__":
    main(sys.argv[1:])
