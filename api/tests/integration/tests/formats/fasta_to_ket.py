﻿import difflib
import os
import sys


def find_diff(a, b):
    return "\n".join(difflib.unified_diff(a.splitlines(), b.splitlines()))


sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa

indigo = Indigo()
indigo.setOption("json-saving-pretty", True)
indigo.setOption("ignore-stereochemistry-errors", True)

print("*** FASTA to KET ***")

root = joinPathPy("molecules/", __file__)
ref_path = joinPathPy("ref/", __file__)

fasta_files = [
    {"file": "test_peptide", "seq_type": "PEPTIDE"},
    {"file": "test_rna", "seq_type": "RNA"},
    {"file": "test_dna", "seq_type": "DNA"},
    {"file": "multiseq", "seq_type": "DNA"},
    {"file": "break", "seq_type": "PEPTIDE"},
    {"file": "comment", "seq_type": "PEPTIDE"},
]

for desc in fasta_files:
    filename = desc["file"]
    mol = indigo.loadFastaFromFile(
        os.path.join(root, filename + ".fasta"), desc["seq_type"]
    )
    # with open(os.path.join(ref_path, filename) + ".ket", "w") as file:
    #     file.write(mol.json())
    with open(os.path.join(ref_path, filename) + ".ket", "r") as file:
        ket_ref = file.read()
    ket = mol.json()
    diff = find_diff(ket_ref, ket)
    if not diff:
        print(filename + ".ket:SUCCEED")
    else:
        print(filename + ".ket:FAILED")
        print(diff)
