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
if canonical.count(".") + 1 != len(SOURCE.split(".")):
    raise Exception("canonical SMILES changed disconnected component multiplicity")
print("component multiplicity: PASS")
print("canonical disconnected components: PASS")

roundtrip = indigo.loadMolecule(canonical)
if roundtrip.canonicalSmiles() != canonical:
    raise Exception("canonical SMILES changed after reload")

permuted_source = ".".join(reversed(SOURCE.split(".")))
permuted = indigo.loadMolecule(permuted_source)
if permuted.canonicalSmiles() != canonical:
    raise Exception("canonical SMILES depends on disconnected component order")
print("component permutation: PASS")
atom_indices = [atom.index() for atom in source.iterateAtoms()]
atom_permuted = source.createSubmolecule(list(reversed(atom_indices)))
if atom_permuted.canonicalSmiles() != canonical:
    raise Exception("canonical SMILES depends on atom order")
print("atom permutation: PASS")

if stereo_count(roundtrip) != source_stereo:
    raise Exception("stereocenter count changed after canonical SMILES roundtrip")
print("stereo preserved: PASS")

components = SOURCE.split(".")
for count in (1, 2, 3, 5, 7):
    selected_components = components[:count]
    repeated_source = ".".join(selected_components)
    reversed_source = ".".join(reversed(selected_components))
    repeated = indigo.loadMolecule(repeated_source).canonicalSmiles()
    reversed_repeated = indigo.loadMolecule(reversed_source).canonicalSmiles()
    if repeated != reversed_repeated or repeated.count(".") + 1 != count:
        raise Exception("canonicalization changed repeated component multiplicity")
print("duplicate multiplicity matrix: PASS")

semantic_pairs = (
    ("stereo", "N[C@@H](C)O.C", "N[C@H](C)O.C"),
    ("isotope", "[13CH4].C", "C.C"),
    ("charge", "[NH4+].O", "N.O"),
    ("bond order", "C=C.O", "CC.O"),
    ("ring size", "C1CC1.O", "C1CCC1.O"),
)
for feature, first, second in semantic_pairs:
    first_smiles = indigo.loadMolecule(first).canonicalSmiles()
    second_smiles = indigo.loadMolecule(second).canonicalSmiles()
    if first_smiles == second_smiles:
        raise Exception("canonical SMILES conflated component " + feature)
print("component semantics: PASS")
aromatic_canonical = indigo.loadMolecule("c1ccccc1.O").canonicalSmiles()
aromatic_reordered = indigo.loadMolecule("O.c1ccccc1").canonicalSmiles()
kekule_canonical = indigo.loadMolecule("C1=CC=CC=C1.O").canonicalSmiles()
kekule_reordered = indigo.loadMolecule("O.C1=CC=CC=C1").canonicalSmiles()
if aromatic_canonical != aromatic_reordered or kekule_canonical != kekule_reordered:
    raise Exception("canonical aromatic components changed after permutation")
if aromatic_canonical == kekule_canonical:
    raise Exception("canonicalization lost aromatic/Kekule state")
if (
    indigo.loadMolecule(aromatic_canonical).canonicalSmiles() != aromatic_canonical
    or indigo.loadMolecule(kekule_canonical).canonicalSmiles() != kekule_canonical
):
    raise Exception("aromatic/Kekule state changed after canonical SMILES roundtrip")
print("aromatic/Kekule state: PASS")

reaction = indigo.loadReaction(
    "[CH3:1][OH:2].[CH3:3][OH:4]>>[CH2:1]=[O:2].[CH2:3]=[O:4]"
)
permuted_reaction = indigo.loadReaction(
    "[CH3:3][OH:4].[CH3:1][OH:2]>>[CH2:3]=[O:4].[CH2:1]=[O:2]"
)
if reaction.canonicalSmiles() != permuted_reaction.canonicalSmiles():
    raise Exception("canonical reaction changed after mapped-component permutation")
print("reaction mappings: PASS")
