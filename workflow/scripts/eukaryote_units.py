"""Reporting units for eukaryotic SSU lineages (SILVA 138.1 rank structure)."""

# SILVA clade names above phylum that carry no identification on their own.
CONTAINERS = {"Animalia", "BCP clade", "Bilateria", "Protostomia", "Deuterostomia", "Lophotrochozoa",
              "Ecdysozoa", "Spiralia", "Ambulacraria", "Panarthropoda", "uncultured", ""}
NEGATIVE_CONTROLS = ("Embryophyta", "Insecta", "Mammalia")
# Ancestors of Metazoa: host reads can resolve this high, so these carry no identification.
METAZOAN_ANCESTORS = ("Amorphea", "Obazoa", "Opisthokonta", "Holozoa", "Choanozoa")


def unit(lineage):
    """Return (kind, reporting unit) for a lineage, or None when it is uninformative.

    kind is one of host, cnidarian, human, negative_control, metazoan, non_metazoan.
    Metazoans are reported by phylum; other eukaryotes by their first five ranks, and
    need at least three resolved ranks.
    """
    ranks = lineage.split(";")
    if ranks[0] != "Eukaryota":
        return None
    if "host siphonophore" in ranks or "Siphonophorae" in ranks:
        return ("host", "Siphonophorae")
    if "Embryophyta" in ranks:
        return ("negative_control", "Embryophyta")
    if "Metazoa" in ranks:
        after = ranks[ranks.index("Metazoa") + 1:]
        informative = [r for r in after if r not in CONTAINERS]
        if not informative:
            return None
        if "Cnidaria" in after:
            return ("cnidarian", "Cnidaria")
        if "Homo" in after:
            return ("human", "Homo")
        for control in ("Insecta", "Mammalia"):
            if control in after:
                return ("negative_control", control)
        return ("metazoan", informative[0])
    if len(ranks) < 4 or ranks[-1] in METAZOAN_ANCESTORS:
        return None
    return ("non_metazoan", ";".join(ranks[1:6]))
