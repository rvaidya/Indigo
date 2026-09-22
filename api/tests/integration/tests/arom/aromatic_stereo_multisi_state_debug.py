from indigo import Indigo


SOURCE = "O1[Si@]2O[Si@@](O[Si@@]3O[Si]1O[Si]O[Si]O3)O[Si]O[Si]O2"
INVALID_AROMATIC = (
    "o1:[Si@]2:o:[Si@@](:o:[Si@@]3:o:[Si]:1:o:[Si]:o:[Si]:o:3)"
    ":o:[Si]:o:[Si]:o:2"
)
ACHIRAL_AROMATIC = INVALID_AROMATIC.replace("@@", "").replace("@", "")


def safe(call):
    try:
        return call()
    except Exception as exc:
        return "ERR: %s" % exc


def dump(label, molecule):
    # Query a clone so diagnostic getters cannot alter the molecule whose state
    # we are trying to describe.
    probe = molecule.clone()
    print("\n=== %s ===" % label)
    print("smiles:", safe(probe.smiles))
    print("stereo count:", len([a for a in probe.iterateStereocenters()]))
    print(
        "idx symbol degree implH valence explicitValence radical charge stereoType"
    )
    for atom in probe.iterateAtoms():
        print(
            atom.atomIndex(),
            safe(atom.symbol),
            safe(atom.degree),
            safe(atom.countImplicitHydrogens),
            safe(atom.valence),
            safe(atom.getExplicitValence),
            safe(atom.radical),
            safe(atom.charge),
            safe(atom.stereocenterType),
        )


indigo = Indigo()

source = indigo.loadMolecule(SOURCE)
dump("SOURCE AFTER LOAD", source)

producer = indigo.loadMolecule(SOURCE)
producer.dearomatize()
dump("SOURCE AFTER DEAROMATIZE", producer)

producer.aromatize()
dump("IN-MEMORY PRODUCER AFTER AROMATIZE", producer)
print("producer canonical:", producer.canonicalSmiles())

fresh = indigo.loadMolecule(ACHIRAL_AROMATIC)
dump("FRESH ACHIRAL LOAD OF SAME AROMATIC GRAPH", fresh)
print("fresh canonical:", fresh.canonicalSmiles())
