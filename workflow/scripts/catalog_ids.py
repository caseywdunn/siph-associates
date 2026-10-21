"""Carry cluster identifiers forward without reusing retired identifiers."""


def stable_ids(
    groups: list[set[str]], previous: dict[str, set[str]], prefix: str, width: int
) -> list[str]:
    """Match clusters by shared members; splits get new IDs and merges keep one.

    Existing groups compete by overlap count, with the old identifier breaking
    ties. A new identifier always exceeds every historical identifier.
    """
    used, result = set(), []
    next_id = max((int(k.removeprefix(prefix)) for k in previous), default=0) + 1
    # Assign the strongest correspondences first so a small split cannot take
    # the old identifier away from the larger surviving part of a cluster.
    matches = sorted(
        (-len(group & members), old, i)
        for i, group in enumerate(groups)
        for old, members in previous.items()
        if group & members
    )
    assigned = {}
    for _, old, i in matches:
        if old not in used and i not in assigned:
            assigned[i] = old
            used.add(old)
    for i in range(len(groups)):
        if i not in assigned:
            assigned[i] = f"{prefix}{next_id:0{width}d}"
            next_id += 1
        result.append(assigned[i])
    if len(result) != len(set(result)):
        raise ValueError("Duplicate catalog identifiers")
    return result
