import os
import sys

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa

indigo = Indigo()

source = "C[C@@H]1[C@@H](S2=N[S@]1=NC(=N2)C(F)(F)F)C"
expected_aromatic = "C[C@H]1[C@@H](C)[s@@]2[n]c([n][s]1[n]2)C(F)(F)F"

mol = indigo.loadMolecule(source)
mol.dearomatize()
mol.aromatize()

aromatic = mol.canonicalSmiles()
assert aromatic == expected_aromatic

roundtrip = indigo.loadMolecule(aromatic)
assert roundtrip.canonicalSmiles() == aromatic

original_stereo = len([atom for atom in mol.iterateStereocenters()])
roundtrip_stereo = len([atom for atom in roundtrip.iterateStereocenters()])
assert roundtrip_stereo == original_stereo

# Benign CX stereo metadata must preserve the explicitly chiral aromatic center.
cx_explicit_group = indigo.loadMolecule(expected_aromatic + " |a:4|")
assert len([atom for atom in cx_explicit_group.iterateStereocenters()]) == original_stereo


# A chemistry edit in another component must trigger final revalidation without
# invalidating an aromatic center whose own chemistry is unchanged. The aromatic
# component has atoms 0-13; the appended star is atom 14 and becomes an R-site.
cx_unrelated_rsite = expected_aromatic + ".[*] |$" + ";" * 14 + "_R1$|"
cx_unrelated = indigo.loadMolecule(cx_unrelated_rsite)
assert len([atom for atom in cx_unrelated.iterateStereocenters()]) == original_stereo


# Multiple explicit aromatic stereocenters must be validated in one load without
# losing either center or depending on component order.
multi_aromatic = expected_aromatic + "." + expected_aromatic
multi = indigo.loadMolecule(multi_aromatic)
multi_saved = multi.canonicalSmiles()
multi_roundtrip = indigo.loadMolecule(multi_saved)
assert len([atom for atom in multi.iterateStereocenters()]) == original_stereo * 2
assert len([atom for atom in multi_roundtrip.iterateStereocenters()]) == original_stereo * 2

connectivity_source = "CC1C(S2=NC(=NS1=N2)C(F)(F)F)C"
connectivity = indigo.loadMolecule(connectivity_source)
connectivity.dearomatize()
connectivity.aromatize()

connectivity_aromatic = connectivity.canonicalSmiles()
connectivity_roundtrip = indigo.loadMolecule(connectivity_aromatic)
assert connectivity_roundtrip.canonicalSmiles() == connectivity_aromatic

invalid_aromatic = "C[c@]1ccccc1"

for invalid in (
    invalid_aromatic,
    invalid_aromatic + "." + expected_aromatic,
    expected_aromatic + "." + invalid_aromatic,
):
    try:
        indigo.loadMolecule(invalid)
    except IndigoException:
        pass
    else:
        raise AssertionError("invalid aromatic carbon chirality was accepted")


cx_group_only = expected_aromatic.replace("[s@@]", "[s]") + " |a:4|"
try:
    indigo.loadMolecule(cx_group_only)
except IndigoException:
    pass
else:
    raise AssertionError(
        "CX stereo group without explicit tetrahedral chirality was accepted"
    )


# CX atom labels are applied after the base SMILES stereocenters are constructed.
# Replacing the validated aromatic sulfur with an R-site must invalidate the
# existing center rather than trusting the pre-CX chemistry.
cx_mutated_center = expected_aromatic + " |$;;;;_R1$|"
try:
    indigo.loadMolecule(cx_mutated_center)
except IndigoException:
    pass
else:
    raise AssertionError("CX-mutated aromatic stereocenter was accepted")


# Tolerant loading must not retain stereo that final CX chemistry invalidates.
indigo.setOption("ignore-stereochemistry-errors", True)
try:
    cx_ignored = indigo.loadMolecule(cx_mutated_center)
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)

ignored_stereo = [atom.index() for atom in cx_ignored.iterateStereocenters()]
assert 4 not in ignored_stereo
assert len(ignored_stereo) == original_stereo - 1
