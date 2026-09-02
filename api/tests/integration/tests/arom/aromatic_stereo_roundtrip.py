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

multi_aromatic = aromatic + "." + aromatic
multi_roundtrip = indigo.loadMolecule(multi_aromatic)
multi_stereo = len([atom for atom in multi_roundtrip.iterateStereocenters()])
assert multi_stereo == 2 * roundtrip_stereo

multi_canonical = multi_roundtrip.canonicalSmiles()
multi_canonical_roundtrip = indigo.loadMolecule(multi_canonical)
assert (
    len([atom for atom in multi_canonical_roundtrip.iterateStereocenters()])
    == multi_stereo
)

connectivity_source = "CC1C(S2=NC(=NS1=N2)C(F)(F)F)C"
connectivity = indigo.loadMolecule(connectivity_source)
connectivity.dearomatize()
connectivity.aromatize()

connectivity_aromatic = connectivity.canonicalSmiles()
connectivity_roundtrip = indigo.loadMolecule(connectivity_aromatic)
assert connectivity_roundtrip.canonicalSmiles() == connectivity_aromatic

try:
    indigo.loadMolecule("C[c@]1ccccc1")
except IndigoException:
    pass
else:
    raise AssertionError("invalid aromatic carbon chirality was accepted")
