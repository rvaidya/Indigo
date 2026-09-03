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
