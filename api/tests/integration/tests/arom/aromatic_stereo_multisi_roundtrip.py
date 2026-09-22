import os
import sys

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *  # noqa


# PubChem CID 135484931 production regression. The concrete source is valid and
# contains three tetrahedral Si stereocenters. Aromatization previously retained
# those centers in memory while producing an aromatic SMILES that a fresh strict
# Indigo instance could not load.
SOURCE = "O1[Si@]2O[Si@@](O[Si@@]3O[Si]1O[Si]O[Si]O3)O[Si]O[Si]O2"
INVALID_AROMATIC = (
    "o1:[Si@]2:o:[Si@@](:o:[Si@@]3:o:[Si]:1:o:[Si]:o:[Si]:o:3)"
    ":o:[Si]:o:[Si]:o:2"
)


def stereo_count(molecule):
    return len([atom for atom in molecule.iterateStereocenters()])


def stereocenter_incident_bond_orders(molecule):
    return [
        [neighbor.bond().bondOrder() for neighbor in atom.iterateNeighbors()]
        for atom in molecule.iterateStereocenters()
    ]


def assert_strict_rejects(smiles):
    try:
        Indigo().loadMolecule(smiles)
    except IndigoException:
        return
    raise AssertionError("invalid aromatic stereocenter was accepted: %s" % smiles)


def invert_parities(smiles):
    marker = "__DOUBLE_AT__"
    return smiles.replace("@@", marker).replace("@", "@@").replace(marker, "@")


# The representation that reached Pocketbook is not independently loadable.
# Keep the consumer boundary strict; the producer must not emit this value.
assert_strict_rejects(INVALID_AROMATIC)

source = Indigo().loadMolecule(SOURCE)
assert stereo_count(source) == 3
assert stereo_count(Indigo().loadMolecule(source.canonicalSmiles())) == 3

# Exercise the producer path used by canonical SMILES generation. The proposed
# aromatic cage is incompatible with the retained stereo when reconstructed
# from serialized state, so this aromatic component must stay concrete.
producer = Indigo().loadMolecule(SOURCE)
producer.dearomatize()
producer.aromatize()
produced = producer.canonicalSmiles()
assert stereo_count(producer) == 3
assert produced != INVALID_AROMATIC
for bond_orders in stereocenter_incident_bond_orders(producer):
    assert 4 not in bond_orders  # BOND_AROMATIC

consumer = Indigo().loadMolecule(produced)
assert stereo_count(consumer) == 3
assert consumer.canonicalSmiles() == produced

# Repeated aromatization must be stable after the incompatible component is
# suppressed.
producer.aromatize()
assert stereo_count(producer) == 3
assert producer.canonicalSmiles() == produced
for bond_orders in stereocenter_incident_bond_orders(producer):
    assert 4 not in bond_orders

# Suppression must be independent of the particular parity spelling.
opposite = Indigo().loadMolecule(invert_parities(SOURCE))
assert stereo_count(opposite) == 3
opposite.dearomatize()
opposite.aromatize()
opposite_canonical = opposite.canonicalSmiles()
assert opposite_canonical != INVALID_AROMATIC
for bond_orders in stereocenter_incident_bond_orders(opposite):
    assert 4 not in bond_orders
assert stereo_count(Indigo().loadMolecule(opposite_canonical)) == 3

# This is a stereo-serialization constraint, not a silicon/oxygen aromaticity
# ban. With stereo removed, the same cage remains free to aromatize.
achiral_source = SOURCE.replace("@@", "").replace("@", "")
achiral = Indigo().loadMolecule(achiral_source)
assert stereo_count(achiral) == 0
achiral.aromatize()
achiral_canonical = achiral.canonicalSmiles()
assert ":" in achiral_canonical
assert stereo_count(Indigo().loadMolecule(achiral_canonical)) == 0
