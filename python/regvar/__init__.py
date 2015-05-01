def removecommongaps(s1, s2):
    """Remove common gap characters between the two sequences.
    Return s1, s2 with these characters removed.
    """
    if len(s1) != len(s2):
        raise ValueError('Sequences must be same length')
    return (
        ''.join(b1 for b1, b2 in zip(s1, s2) if b1 != '-' or b2 != '-'),
        ''.join(b2 for b1, b2 in zip(s1, s2) if b1 != '-' or b2 != '-'),
    )


def seqdist(s1, s2, mismatchpen=-.5, gapopenpen=-.25, gapextendpen=-.05):
    """
    The distance between two sequences.
    """
    # s1, s2 = removecommongaps(s1, s2)
    from Bio.pairwise2 import align, format_alignment
    alignment = next(iter(align.globalms(
        s1, s2, 1, mismatchpen, gapopenpen, gapextendpen)))
    print(format_alignment(*alignment))
    return alignment[2]


def index_of_nth_base(gappedseq, n):
    """Return the index of the nth non-gapped base in the sequence
    where n=0 is the first."""
    nongapped = 0
    for i, b in enumerate(gappedseq):
        if b == '-':
            continue
        if nongapped == n:
            return i
        nongapped += 1
    raise ValueError(
        "Could not find {0}'th non-gapped base in sequence".format(n))
