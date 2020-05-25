"""
Regression tests for codonharmonizer

"""

import io
import os.path
import subprocess

import pandas as pd
import Bio.SeqIO

def path(x):
    dirn = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(dirn, x)

def test_reference():
    output = subprocess.check_output(["codonharmonizer",
        "--write_freqs",
        path("GCA_000011385.1_ASM1138v1_cds_from_genomic.fna"),
    ])
    output = io.StringIO(output.decode('utf-8'))
    df_test = pd.read_csv(output)
    df_test.sort_values('codon', inplace=True)

    df_expected = pd.read_csv(path("Gvio_Freq.csv"))
    df_expected.sort_values('codon', inplace=True)

    assert df_test.equals(df_expected)
    return

def test_recode():
    output = subprocess.check_output(["codonharmonizer",
        "--source", path("Gvio_Freq.csv"),
        "--target", path("Eco_MG1655_Freq.csv"),
        path("example_gene.faa")
    ])
    output = io.StringIO(output.decode('utf-8'))
    test_recode = Bio.SeqIO.read(output, "fasta")

    expected_recode = Bio.SeqIO.read(path("example_gene.recode.faa"), "fasta")
    assert expected_recode.seq == test_recode.seq
    return
