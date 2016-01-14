def get_indel(aligned_segment, query_position):
    """
    Similar behavior as the indel attribute in the pysam.PileupRead class (0-based query_position)
    MATCH       - 0
    DELETION    - Negative integer on the query position before the deletion for the number of bases
    INSERTION   - 1 on the insertion site
    """
    sum = 0
    for tuple in aligned_segment.cigartuples:
        # (0, 1)
        if tuple[0] != 2:
            sum = sum + tuple[1]
        if query_position < sum:
            pass
    pass
