"""
Regression tests for codonharmonizer

"""

import io
import os.path
import subprocess
import re

import numpy as np
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

    df_expected = pd.read_csv(path("Gvio.freq.csv"))
    df_expected.sort_values('codon', inplace=True)

    assert df_test.equals(df_expected)
    return

def test_recode():
    output = subprocess.check_output(["codonharmonizer",
        "--source", path("Gvio.freq.csv"),
        "--target", path("Ecol_MG1655.freq.csv"),
        path("example_gene.fna")
    ])
    output = io.StringIO(output.decode('utf-8'))
    test_recode = Bio.SeqIO.read(output, "fasta")

    # Check harmonized sequence
    expected_recode = Bio.SeqIO.read(path("example_gene.recode.fna"), "fasta")
    assert expected_recode.seq == test_recode.seq

    # Check calculated CHI values
    expected_initial_chi = 0.2219
    expected_harmonized_chi = 0.0991

    # last bit of header
    m = re.search(r'CHI_0:(\d+\.\d+).*CHI:(\d+\.\d+)', test_recode.description)
    test_initial_chi = float(m.group(1))
    test_harmonized_chi = float(m.group(2))

    assert np.isclose(expected_initial_chi, test_initial_chi)
    assert np.isclose(expected_harmonized_chi, test_harmonized_chi)

    return
