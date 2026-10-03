import os
import sys

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa


SOURCE = (
    "CC(C)(C)[C@@H]1CCCCC[C@H]1C(C)(C)C.CC(C)(C)[C@@H]1CCCCC[C@H]1C(C)(C)C."
    "CC(C)(C)[C@@H]1CCCCC[C@@H]1C(C)(C)C.CC(C)(C)[C@H]1CCCCC[C@@H]1C(C)(C)C."
    "CC(C)(C)[C@H]1CCCCC[C@@H]1C(C)(C)C.CC(C)(C)[C@@H]1CCCC[C@H]1C(C)(C)C."
    "CC(C)(C)[C@@H]1CCCC[C@H]1C(C)(C)C.CC(C)(C)[C@@H]1CCCC[C@H]1C(C)(C)C."
    "CC(C)(C)[C@@H]1CCCC[C@@H]1C(C)(C)C.CC(C)(C)[C@H]1CCCC[C@@H]1C(C)(C)C."
    "CC(C)(C)[C@H]1CCCC[C@@H]1C(C)(C)C.CC(C)(C)[C@H]1CCCC[C@@H]1C(C)(C)C."
    "CC(C)(C)[C@@H]1CCC[C@H]1C(C)(C)C.CC(C)(C)[C@@H]1CCC[C@H]1C(C)(C)C."
    "CC(C)(C)[C@@H]1CCC[C@H]1C(C)(C)C.CC(C)(C)[C@@H]1CCC[C@@H]1C(C)(C)C."
    "CC(C)(C)[C@H]1CCC[C@@H]1C(C)(C)C.CC(C)(C)[C@H]1CCC[C@@H]1C(C)(C)C."
    "CC(C)(C)[C@H]1CCC[C@@H]1C(C)(C)C"
)


def stereo_count(mol):
    return len(list(mol.iterateStereocenters()))


indigo = Indigo()
indigo.setOption("timeout", "5000")

source = indigo.loadMolecule(SOURCE)
source_stereo = stereo_count(source)
canonical = source.canonicalSmiles()
print("canonical disconnected components: PASS")

roundtrip = indigo.loadMolecule(canonical)
if roundtrip.canonicalSmiles() != canonical:
    raise Exception("canonical SMILES changed after reload")
print("roundtrip: PASS")

permuted_source = ".".join(reversed(SOURCE.split(".")))
permuted = indigo.loadMolecule(permuted_source)
if permuted.canonicalSmiles() != canonical:
    raise Exception("canonical SMILES depends on disconnected component order")
print("component permutation: PASS")

if stereo_count(roundtrip) != source_stereo:
    raise Exception("stereocenter count changed after canonical SMILES roundtrip")
print("stereo preserved: PASS")
