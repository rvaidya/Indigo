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

# The opposite sulfur parity must remain a distinct stereoisomer and round-trip.
opposite_aromatic = expected_aromatic.replace("[s@@]", "[s@]", 1)
opposite = indigo.loadMolecule(opposite_aromatic)
opposite_saved = opposite.canonicalSmiles()
assert opposite_saved != aromatic
assert len([atom for atom in opposite.iterateStereocenters()]) == original_stereo
assert indigo.loadMolecule(opposite_saved).canonicalSmiles() == opposite_saved

# Two explicit aromatic sulfur centers in the same aromatic system exercise
# joint Kekule compatibility rather than independent component validation.
joint_source = "C[C@@H]1[C@@H]([S@]2=N[S@]1=NC(=N2)C(F)(F)F)C"
joint = indigo.loadMolecule(joint_source)
joint.dearomatize()
joint.aromatize()
joint_aromatic = joint.canonicalSmiles()
joint_stereo = len([atom for atom in joint.iterateStereocenters()])
joint_roundtrip = indigo.loadMolecule(joint_aromatic)
assert joint_stereo == original_stereo + 1
assert len([atom for atom in joint_roundtrip.iterateStereocenters()]) == joint_stereo
assert joint_roundtrip.canonicalSmiles() == joint_aromatic

# Benign CX stereo metadata must preserve the explicitly chiral aromatic center.
cx_explicit_group = indigo.loadMolecule(expected_aromatic + " |a:4|")
assert len([atom for atom in cx_explicit_group.iterateStereocenters()]) == original_stereo


# A chemistry edit in another component must trigger final revalidation without
# invalidating an aromatic center whose own chemistry is unchanged. The aromatic
# component has atoms 0-13; the appended star is atom 14 and becomes an R-site.
cx_unrelated_rsite = expected_aromatic + ".[*] |$" + ";" * 14 + "_R1$|"
cx_unrelated = indigo.loadMolecule(cx_unrelated_rsite)
assert len([atom for atom in cx_unrelated.iterateStereocenters()]) == original_stereo

# Coordinates advance the molecule edit revision and rebuild bond-derived stereo.
# This must trigger final aromatic revalidation without changing the chemistry.
coordinates = ";".join("%d,%d," % (i, i % 3) for i in range(14))
cx_coordinates = indigo.loadMolecule(expected_aromatic + " |(" + coordinates + ")|")
assert cx_coordinates.canonicalSmiles() == aromatic
assert len([atom for atom in cx_coordinates.iterateStereocenters()]) == original_stereo


# Multiple explicit aromatic stereocenters must be validated in one load without
# losing either center or depending on component order.
multi_aromatic = expected_aromatic + "." + expected_aromatic
multi = indigo.loadMolecule(multi_aromatic)
multi_saved = multi.canonicalSmiles()
multi_roundtrip = indigo.loadMolecule(multi_saved)
assert len([atom for atom in multi.iterateStereocenters()]) == original_stereo * 2
assert len([atom for atom in multi_roundtrip.iterateStereocenters()]) == original_stereo * 2

# Finite CurlySMILES repetition copies stereocenters through a submolecule
# mapping. Aromatic fallback provenance must follow those copied centers so the
# final validator checks them instead of rejecting or silently bypassing them.
curly_repeated = expected_aromatic.replace("C[C@H]1", "C{-}[C@H]1", 1) + "{+nn=2}"
curly = indigo.loadMolecule(curly_repeated)
assert len([atom for atom in curly.iterateStereocenters()]) == original_stereo * 2

curly_repeated_three = curly_repeated.replace("{+nn=2}", "{+nn=3}")
curly_three = indigo.loadMolecule(curly_repeated_three)
assert len([atom for atom in curly_three.iterateStereocenters()]) == original_stereo * 3

indigo.setOption("ignore-stereochemistry-errors", True)
try:
    curly_tolerant = indigo.loadMolecule(curly_repeated)
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)

assert len([atom for atom in curly_tolerant.iterateStereocenters()]) == original_stereo * 2

connectivity_source = "CC1C(S2=NC(=NS1=N2)C(F)(F)F)C"
connectivity = indigo.loadMolecule(connectivity_source)
connectivity.dearomatize()
connectivity.aromatize()

connectivity_aromatic = connectivity.canonicalSmiles()
connectivity_roundtrip = indigo.loadMolecule(connectivity_aromatic)
assert connectivity_roundtrip.canonicalSmiles() == connectivity_aromatic

invalid_aromatic = "C[c@]1ccccc1"

# The fallback is sulfur-specific. A normal N-substituted pyrrole-like
# aromatic system is valid, but explicit tetrahedral chirality on that
# aromatic nitrogen must not become legal through Kekule fallback.
neutral_aromatic_n = "Cn1cccc1"
indigo.loadMolecule(neutral_aromatic_n)

invalid_aromatic_n = ("C[n@]1cccc1", "C[n@@]1cccc1")
for invalid in invalid_aromatic_n:
    try:
        indigo.loadMolecule(invalid)
    except IndigoException:
        pass
    else:
        raise AssertionError("aromatic nitrogen chirality was accepted")

indigo.setOption("ignore-stereochemistry-errors", True)
try:
    tolerant_aromatic_n = [indigo.loadMolecule(value) for value in invalid_aromatic_n]
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)

for molecule in tolerant_aromatic_n:
    assert len([atom for atom in molecule.iterateStereocenters()]) == 0

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


# Tolerant mode must keep valid aromatic fallback stereo while omitting
# aromatic chirality that has no globally valid Kekule configuration.
indigo.setOption("ignore-stereochemistry-errors", True)
try:
    tolerant_valid = indigo.loadMolecule(expected_aromatic)
    tolerant_invalid = indigo.loadMolecule(invalid_aromatic)
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)

assert len([atom for atom in tolerant_valid.iterateStereocenters()]) == original_stereo
assert len([atom for atom in tolerant_invalid.iterateStereocenters()]) == 0


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
