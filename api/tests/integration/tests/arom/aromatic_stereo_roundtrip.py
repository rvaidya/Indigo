import os
import sys

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa

indigo = Indigo()


def stereo_count(molecule):
    return len([atom for atom in molecule.iterateStereocenters()])


def assert_strict_canonical_roundtrip(smiles, expected_stereo=None):
    # Use a fresh Indigo session for the consumer side. The contract being
    # tested is serializer -> independent loader, not merely repeated calls on
    # one molecule object.
    producer = Indigo()
    molecule = producer.loadMolecule(smiles)
    saved = molecule.canonicalSmiles()
    if expected_stereo is not None:
        assert stereo_count(molecule) == expected_stereo

    consumer = Indigo()
    reloaded = consumer.loadMolecule(saved)
    if expected_stereo is not None:
        assert stereo_count(reloaded) == expected_stereo
    assert reloaded.canonicalSmiles() == saved
    return saved


def assert_aromatic_serializer_contract(source_smiles, aromatic_atom_token):
    # Start from an explicit/Kekule representation that the ordinary Indigo
    # stereocenter rules already accept. If Indigo aromatizes and serializes
    # that center as lowercase @/@@, a fresh Indigo loader must consume the
    # emitted representation without losing stereochemistry.
    producer = Indigo()
    molecule = producer.loadMolecule(source_smiles)
    source_canonical = molecule.canonicalSmiles()
    expected_stereo = stereo_count(molecule)
    assert expected_stereo > 0

    molecule.dearomatize()
    assert stereo_count(molecule) == expected_stereo
    molecule.aromatize()
    assert stereo_count(molecule) == expected_stereo

    aromatic = molecule.canonicalSmiles()
    assert aromatic_atom_token in aromatic

    consumer = Indigo()
    reloaded = consumer.loadMolecule(aromatic)
    assert stereo_count(reloaded) == expected_stereo
    assert reloaded.canonicalSmiles() == aromatic

    reloaded.dearomatize()
    assert stereo_count(reloaded) == expected_stereo
    assert reloaded.canonicalSmiles() == source_canonical
    return aromatic, expected_stereo


source = "C[C@@H]1[C@@H](S2=N[S@]1=NC(=N2)C(F)(F)F)C"
expected_aromatic = "C[C@H]1[C@@H](C)[s@@]2[n]c([n][s]1[n]2)C(F)(F)F"

mol = indigo.loadMolecule(source)
source_canonical = mol.canonicalSmiles()
original_stereo = len([atom for atom in mol.iterateStereocenters()])

mol.dearomatize()
assert len([atom for atom in mol.iterateStereocenters()]) == original_stereo
mol.aromatize()
assert len([atom for atom in mol.iterateStereocenters()]) == original_stereo

aromatic = mol.canonicalSmiles()
assert aromatic == expected_aromatic

roundtrip = indigo.loadMolecule(aromatic)
assert roundtrip.canonicalSmiles() == aromatic

roundtrip_stereo = len([atom for atom in roundtrip.iterateStereocenters()])
assert roundtrip_stereo == original_stereo

roundtrip_dearomatized = indigo.loadMolecule(aromatic)
roundtrip_dearomatized.dearomatize()
assert len([atom for atom in roundtrip_dearomatized.iterateStereocenters()]) == original_stereo
assert roundtrip_dearomatized.canonicalSmiles() == source_canonical

# The regression is the serializer/loader contract, not a sulfur-only special
# case. These explicit/Kekule sources exercise the heavy aromatic elements for
# which Indigo already defines tetrahedral stereocenter chemistry and can emit
# lowercase aromatic SMILES.
serializer_contract_cases = (
    (
        source,
        "[s@",
    ),
    (
        "CN(C)[P@]1(F)=NP(F)(F)=NP(=N1)(N(C)C)N(C)C",
        "[p@",
    ),
    (
        "C[As@]1(F)C=CC=C1",
        "[as@",
    ),
)
serializer_contract_outputs = []
for contract_source, aromatic_token in serializer_contract_cases:
    serializer_contract_outputs.append(
        assert_aromatic_serializer_contract(contract_source, aromatic_token)
    )

# CID 16419269 exposed the phosphorus member of the same bug class after the
# sulfur fix had already passed the full release harness. Keep the exact
# Indigo-generated aromatic canonical SMILES as a permanent production
# regression in addition to the producer-side contract above.
pubchem_16419269_aromatic = (
    "CN(C)[p@]1(F)[n][p](F)(F)[n][p]([n]1)(N(C)C)N(C)C"
)
pubchem_16419269 = Indigo().loadMolecule(pubchem_16419269_aromatic)
assert stereo_count(pubchem_16419269) == 1
pubchem_16419269_saved = pubchem_16419269.canonicalSmiles()
assert_strict_canonical_roundtrip(pubchem_16419269_saved, 1)

# The opposite phosphorus parity is a distinct stereoisomer and must satisfy
# the same independent save/reload contract.
pubchem_16419269_opposite = pubchem_16419269_aromatic.replace("[p@]", "[p@@]", 1)
pubchem_16419269_opposite_saved = assert_strict_canonical_roundtrip(
    pubchem_16419269_opposite, 1
)
assert pubchem_16419269_opposite_saved != pubchem_16419269_saved

# Different supported heavy-element centers must coexist in one molecule load;
# validation context must not accidentally depend on all fallback centers being
# the same element.
mixed_heavy_aromatic = expected_aromatic + "." + pubchem_16419269_aromatic
mixed_heavy = Indigo().loadMolecule(mixed_heavy_aromatic)
mixed_heavy_stereo = original_stereo + 1
assert stereo_count(mixed_heavy) == mixed_heavy_stereo
mixed_heavy_saved = mixed_heavy.canonicalSmiles()
assert_strict_canonical_roundtrip(mixed_heavy_saved, mixed_heavy_stereo)

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
assert multi_roundtrip.canonicalSmiles() == multi_saved

# Finite CurlySMILES repetition copies stereocenters through a submolecule
# mapping. Aromatic fallback provenance must follow those copied centers so the
# final validator checks them instead of rejecting or silently bypassing them.
curly_repeated = expected_aromatic.replace("C[C@H]1", "C{-}[C@H]1", 1) + "{+nn=2}"
curly = indigo.loadMolecule(curly_repeated)
assert len([atom for atom in curly.iterateStereocenters()]) == original_stereo * 2
curly_saved = curly.canonicalSmiles()
curly_roundtrip = indigo.loadMolecule(curly_saved)
assert len([atom for atom in curly_roundtrip.iterateStereocenters()]) == original_stereo * 2
assert curly_roundtrip.canonicalSmiles() == curly_saved

curly_repeated_three = curly_repeated.replace("{+nn=2}", "{+nn=3}")
curly_three = indigo.loadMolecule(curly_repeated_three)
assert len([atom for atom in curly_three.iterateStereocenters()]) == original_stereo * 3
curly_three_saved = curly_three.canonicalSmiles()
curly_three_roundtrip = indigo.loadMolecule(curly_three_saved)
assert len([atom for atom in curly_three_roundtrip.iterateStereocenters()]) == original_stereo * 3
assert curly_three_roundtrip.canonicalSmiles() == curly_three_saved

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

# The fallback is deliberately limited to the heavy-element stereocenter
# family above. A normal N-substituted pyrrole-like aromatic system is valid,
# but explicit tetrahedral chirality on that aromatic nitrogen must not become
# legal through Kekule fallback.
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
    assert stereo_count(molecule) == 0

# Selenium is lowercase-aromatic-capable in the SMILES saver, but Indigo has no
# corresponding tetrahedral stereocenter configuration. The compatibility
# fallback must therefore not turn aromatic selenium into a new stereo class.
invalid_aromatic_se = ("C[se@]1cccc1", "C[se@@]1cccc1")
for invalid in invalid_aromatic_se:
    try:
        Indigo().loadMolecule(invalid)
    except IndigoException:
        pass
    else:
        raise AssertionError("aromatic selenium chirality was accepted")

# An unrelated neutral aromatic nitrogen component must not interfere with the
# valid heavy-element fallback, regardless of component order.
for combined in (
    expected_aromatic + "." + neutral_aromatic_n,
    neutral_aromatic_n + "." + expected_aromatic,
):
    combined_molecule = indigo.loadMolecule(combined)
    combined_saved = combined_molecule.canonicalSmiles()
    assert len([atom for atom in combined_molecule.iterateStereocenters()]) == original_stereo
    combined_roundtrip = indigo.loadMolecule(combined_saved)
    assert len([atom for atom in combined_roundtrip.iterateStereocenters()]) == original_stereo
    assert combined_roundtrip.canonicalSmiles() == combined_saved

# In tolerant mode an invalid aromatic N center must disappear while the valid
# supported heavy-element stereo survives, in either component order. The sanitized save must
# then reload strictly with the same stereo state.
mixed_tolerant_inputs = (
    expected_aromatic + "." + invalid_aromatic_n[0],
    invalid_aromatic_n[0] + "." + expected_aromatic,
)
indigo.setOption("ignore-stereochemistry-errors", True)
try:
    mixed_tolerant = [indigo.loadMolecule(value) for value in mixed_tolerant_inputs]
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)

for molecule in mixed_tolerant:
    assert len([atom for atom in molecule.iterateStereocenters()]) == original_stereo
    mixed_saved = molecule.canonicalSmiles()
    mixed_roundtrip = indigo.loadMolecule(mixed_saved)
    assert len([atom for atom in mixed_roundtrip.iterateStereocenters()]) == original_stereo
    assert mixed_roundtrip.canonicalSmiles() == mixed_saved

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
