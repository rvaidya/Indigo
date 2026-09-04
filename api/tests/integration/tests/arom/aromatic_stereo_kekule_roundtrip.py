import os
import sys

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa


PUBCHEM_S_SOURCE = "C[C@@H]1[C@@H](S2=N[S@]1=NC(=N2)C(F)(F)F)C"
PUBCHEM_S_AROMATIC = "C[C@H]1[C@@H](C)[s@@]2[n]c([n][s]1[n]2)C(F)(F)F"

# Current PubChem source connectivity for CID 21732325. It intentionally has
# no defined stereochemistry at phosphorus.
PUBCHEM_P_SOURCE = "CN(C)P1(=NP(=NP(=N1)(F)F)(N(C)C)F)N(C)C"

# Exact stored Indigo aromatic representation that failed to reload during
# PubChem reconciliation for this compound (internal compounds.id 16419269).
PUBCHEM_P_AROMATIC = "CN(C)[p@]1(F)[n][p](F)(F)[n][p]([n]1)(N(C)C)N(C)C"


def stereo_count(molecule):
    return len([atom for atom in molecule.iterateStereocenters()])


def assert_canonical_roundtrip(smiles, expected_stereo=None):
    producer = Indigo()
    molecule = producer.loadMolecule(smiles)
    canonical = molecule.canonicalSmiles()

    if expected_stereo is not None:
        assert stereo_count(molecule) == expected_stereo
    assert "@" in canonical

    consumer = Indigo()
    reloaded = consumer.loadMolecule(canonical)
    if expected_stereo is not None:
        assert stereo_count(reloaded) == expected_stereo
    assert reloaded.canonicalSmiles() == canonical
    return canonical


def assert_aromatic_serializer_contract(source, aromatic_marker):
    producer = Indigo()
    molecule = producer.loadMolecule(source)
    source_stereo = stereo_count(molecule)
    assert source_stereo > 0

    # The bug class is created by retaining a valid ordinary stereocenter while
    # its ring bonds become aromatic. Exercise the real producer path.
    molecule.dearomatize()
    molecule.aromatize()
    aromatic = molecule.canonicalSmiles()

    assert aromatic_marker in aromatic
    assert "@" in aromatic
    assert stereo_count(molecule) == source_stereo

    consumer = Indigo()
    reloaded = consumer.loadMolecule(aromatic)
    assert stereo_count(reloaded) == source_stereo
    assert reloaded.canonicalSmiles() == aromatic
    return aromatic, source_stereo


def assert_strict_rejects(smiles):
    try:
        Indigo().loadMolecule(smiles)
    except IndigoException:
        return
    raise AssertionError("invalid aromatic stereocenter was accepted: %s" % smiles)


def assert_tolerant_omits_stereo(smiles):
    indigo = Indigo()
    indigo.setOption("ignore-stereochemistry-errors", True)
    try:
        molecule = indigo.loadMolecule(smiles)
    finally:
        indigo.setOption("ignore-stereochemistry-errors", False)

    assert stereo_count(molecule) == 0
    return molecule


# ---------------------------------------------------------------------------
# Exact production regressions
# ---------------------------------------------------------------------------

s_aromatic, s_stereo = assert_aromatic_serializer_contract(
    PUBCHEM_S_SOURCE, "[s@"
)
assert s_aromatic == PUBCHEM_S_AROMATIC

s_dearomatized = Indigo().loadMolecule(s_aromatic)
s_dearomatized.dearomatize()
s_source_canonical = Indigo().loadMolecule(PUBCHEM_S_SOURCE).canonicalSmiles()
assert s_dearomatized.canonicalSmiles() == s_source_canonical

# CID 21732325's current PubChem source does not define P stereo, so it is not
# a parity baseline. It is still an independent connectivity baseline: Indigo
# must not invent stereo while canonicalizing it, and the historical aromatic
# failure must reduce to the same source chemistry when its explicit stereo is
# removed.
p_source = Indigo().loadMolecule(PUBCHEM_P_SOURCE)
assert stereo_count(p_source) == 0
p_source.dearomatize()
p_source_canonical = p_source.canonicalSmiles()
assert "@" not in p_source_canonical

p_source.aromatize()
p_source_aromatic = p_source.canonicalSmiles()
assert stereo_count(p_source) == 0
assert "@" not in p_source_aromatic

p_molecule = Indigo().loadMolecule(PUBCHEM_P_AROMATIC)
p_stereo = stereo_count(p_molecule)
assert p_stereo == 1
p_aromatic = p_molecule.canonicalSmiles()
assert p_aromatic == PUBCHEM_P_AROMATIC
assert_canonical_roundtrip(PUBCHEM_P_AROMATIC, p_stereo)

# Both parities must remain distinct and independently stable.
s_opposite = PUBCHEM_S_AROMATIC.replace("[s@@]", "[s@]", 1)
p_opposite = PUBCHEM_P_AROMATIC.replace("[p@]", "[p@@]", 1)
assert assert_canonical_roundtrip(s_opposite, s_stereo) != s_aromatic
assert assert_canonical_roundtrip(p_opposite, p_stereo) != p_aromatic

# The source does not choose between these P enantiomeric representations.
# Removing explicit stereo from either must recover the same PubChem source
# connectivity/bonding.
for p_variant in (PUBCHEM_P_AROMATIC, p_opposite):
    p_variant_without_stereo = Indigo().loadMolecule(p_variant)
    p_variant_without_stereo.clearStereocenters()
    assert p_variant_without_stereo.canonicalSmiles() == p_source_aromatic

    p_variant_without_stereo.dearomatize()
    assert p_variant_without_stereo.canonicalSmiles() == p_source_canonical


# ---------------------------------------------------------------------------
# Stereocenter configuration classes that aromatization can expose
# ---------------------------------------------------------------------------

# These are not element whitelists. Each source is accepted by Indigo's normal
# Kekule stereocenter rules first, then aromatized. They exercise zero-double
# pyramidal centers where the aromatic ring can use a lone-pair state.
zero_double_sources = (
    ("neutral nitrogen", "C[N@]1C=C(C)C=C1", "[n@"),
    ("neutral phosphorus", "C[P@]1C=C(C)C=C1", "[p@"),
    ("sulfonium sulfur", "C[S@+]1C=C(C)C=C1", "[s@"),
)
for _, source, marker in zero_double_sources:
    assert_aromatic_serializer_contract(source, marker)

# One-double configurations are a separate stereocenter-table class. The real
# S/P regressions above cover neutral S and neutral P; these additional cases
# exercise neutral pentavalent N and tetravalent S+ without adding loader
# element policy. Use a six-member alternating ring so every ring atom can
# participate in the aromatic system; the methyl branch also makes the two ring
# paths around the stereocenter distinct.
one_double_sources = (
    ("neutral pentavalent nitrogen", "C[N@]1(F)=C(C)C=CC=C1", "[n@"),
    ("tetravalent sulfonium sulfur", "C[S@+]1(F)=C(C)C=CC=C1", "[s@"),
)
for _, source, marker in one_double_sources:
    assert_aromatic_serializer_contract(source, marker)


# ---------------------------------------------------------------------------
# Boundary cases: local stereocenter rules are not enough
# ---------------------------------------------------------------------------

# Aromatic carbon can look tetrahedral only if both incident aromatic bonds are
# made single locally. No globally valid benzene-like Kekule state permits that
# assignment, so it must still be rejected.
invalid_aromatic_carbon = "C[c@]1cc(C)cc1"
assert_strict_rejects(invalid_aromatic_carbon)
assert_tolerant_omits_stereo(invalid_aromatic_carbon)

# In tolerant mode an invalid aromatic center must be omitted without removing
# a valid fallback center in another component. The sanitized result must then
# be consumable strictly.
for mixed_tolerant_input in (
    PUBCHEM_S_AROMATIC + "." + invalid_aromatic_carbon,
    invalid_aromatic_carbon + "." + PUBCHEM_S_AROMATIC,
):
    indigo = Indigo()
    indigo.setOption("ignore-stereochemistry-errors", True)
    try:
        mixed_tolerant = indigo.loadMolecule(mixed_tolerant_input)
    finally:
        indigo.setOption("ignore-stereochemistry-errors", False)

    assert stereo_count(mixed_tolerant) == s_stereo
    mixed_tolerant_canonical = mixed_tolerant.canonicalSmiles()

    # Tolerant cleanup must produce a structure that is independently valid in
    # strict mode, regardless of which component was encountered first.
    strict_mixed_tolerant = Indigo().loadMolecule(mixed_tolerant_canonical)
    assert stereo_count(strict_mixed_tolerant) == s_stereo

# Indigo already has aromatic silicon fixtures whose canonical form uses
# uppercase Si plus explicit aromatic ':' bonds. A Si tetrahedral configuration
# requires zero double bonds, but this aromatic ring requires Si to participate
# in the Kekule matching. This specifically proves that fallback detection uses
# actual aromatic bonds rather than lowercase atom spelling.
invalid_aromatic_silicon = "C[Si@]1:ccccc:1"
assert_strict_rejects(invalid_aromatic_silicon)
assert_tolerant_omits_stereo(invalid_aromatic_silicon)

# Tetracoordinate B- is an ordinary Indigo stereocenter class, but forcing it
# into this aromatic ring does not by itself make a globally valid aromatic
# stereocenter. This guards against replacing an element whitelist with
# unconditional "any aromatic atom" acceptance.
invalid_aromatic_boron = "C[B@-]1(F):c(C):c:c:c:1"
assert_strict_rejects(invalid_aromatic_boron)

# Se and Te are aromatic-capable in Indigo, but Indigo's current tetrahedral
# stereocenter table has no matching configurations. The generic fallback must
# therefore fail through the existing stereocenter model, not a loader list.
for invalid in (
    "C[se@]1cc(C)cc1",
    "C[te@]1cc(C)cc1",
):
    assert_strict_rejects(invalid)
    assert_tolerant_omits_stereo(invalid)


# ---------------------------------------------------------------------------
# Multiple centers and global compatibility
# ---------------------------------------------------------------------------

# Different configuration classes/elements must coexist without sharing hidden
# element-specific state.
mixed_stereo = s_stereo + p_stereo
for mixed_aromatic in (
    PUBCHEM_S_AROMATIC + "." + PUBCHEM_P_AROMATIC,
    PUBCHEM_P_AROMATIC + "." + PUBCHEM_S_AROMATIC,
):
    mixed = Indigo().loadMolecule(mixed_aromatic)
    assert stereo_count(mixed) == mixed_stereo
    assert_canonical_roundtrip(mixed.canonicalSmiles(), mixed_stereo)

# Two explicit aromatic centers in the same ring system must be satisfiable in
# one Kekule assignment, not validated independently.
joint_source = "C[C@@H]1[C@@H]([S@]2=N[S@]1=NC(=N2)C(F)(F)F)C"
joint = Indigo().loadMolecule(joint_source)
joint.dearomatize()
joint.aromatize()
joint_aromatic = joint.canonicalSmiles()
joint_stereo = stereo_count(joint)
assert joint_stereo == s_stereo + 1
assert_canonical_roundtrip(joint_aromatic, joint_stereo)


# ---------------------------------------------------------------------------
# Post-SMILES chemistry changes must not leave stale fallback stereo
# ---------------------------------------------------------------------------

# Benign coordinates edit molecule state but not chemistry; stereo must survive.
coordinates = ";".join("%d,%d," % (i, i % 3) for i in range(14))
with_coordinates = Indigo().loadMolecule(
    PUBCHEM_S_AROMATIC + " |(" + coordinates + ")|"
)
assert stereo_count(with_coordinates) == s_stereo
assert with_coordinates.canonicalSmiles() == s_aromatic

# Replacing the validated sulfur with an R-site after base stereo construction
# invalidates the center. Strict mode rejects; tolerant mode removes only the
# stale fallback center.
cx_mutated_s = PUBCHEM_S_AROMATIC + " |$;;;;_R1$|"
assert_strict_rejects(cx_mutated_s)

indigo = Indigo()
indigo.setOption("ignore-stereochemistry-errors", True)
try:
    cx_tolerant_s = indigo.loadMolecule(cx_mutated_s)
finally:
    indigo.setOption("ignore-stereochemistry-errors", False)
assert 4 not in [atom.index() for atom in cx_tolerant_s.iterateStereocenters()]
assert stereo_count(cx_tolerant_s) == s_stereo - 1

# ---------------------------------------------------------------------------
# CurlySMILES copies must retain fallback provenance and revalidate every copy
# ---------------------------------------------------------------------------

curly_two = PUBCHEM_S_AROMATIC.replace(
    "C[C@H]1", "C{-}[C@H]1", 1
) + "{+nn=2}"
curly_two_molecule = Indigo().loadMolecule(curly_two)
assert stereo_count(curly_two_molecule) == s_stereo * 2
assert_canonical_roundtrip(
    curly_two_molecule.canonicalSmiles(), s_stereo * 2
)

curly_three = curly_two.replace("{+nn=2}", "{+nn=3}")
curly_three_molecule = Indigo().loadMolecule(curly_three)
assert stereo_count(curly_three_molecule) == s_stereo * 3
assert_canonical_roundtrip(
    curly_three_molecule.canonicalSmiles(), s_stereo * 3
)


# Query/SMARTS remain source-isolated because the exceptional validator requires
# a concrete Molecule (_mol != nullptr); no query-specific behavior is changed.

# Ordinary aromatic and ordinary tetrahedral loading remain stable.
assert Indigo().loadMolecule("Cn1cccc1").canonicalSmiles()
assert_canonical_roundtrip("N[C@@H](C)C(=O)O", 1)
