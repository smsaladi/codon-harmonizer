"""
Regression tests for codonharmonizer

"""

import io
import subprocess

import pandas as pd
import Bio.SeqIO

def test_reference():
    output = subprocess.check_output(["codonharmonizer",
        "--write_freqs",
        "GCA_000011385.1_ASM1138v1_cds_from_genomic.fna",
    ])
    output = io.StringIO(output.decode('utf-8'))
    df_test = pd.read_csv(output)
    df_test.sort_values('codon', inplace=True)

    df_expected = pd.read_csv("Gvio_Freq.csv")
    df_expected.sort_values('codon', inplace=True)

    assert df_test.equals(df_expected)
    return

def test_recode():
    output = subprocess.check_output(["codonharmonizer",
        "--source", "Gvio_Freq.csv",
        "--target", "Eco_MG1655_Freq.csv",
        "example_gene.faa"
    ])
    output = io.StringIO(output.decode('utf-8'))
    test_recode = Bio.SeqIO.read(output, "fasta")

    expected_recode = Bio.SeqIO.read("example_gene.recode.faa", "fasta")
    assert expected_recode.seq == test_recode.seq
    return
