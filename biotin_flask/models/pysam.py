def get_indel(aligned_segment, query_position):
    """
    Similar behavior as the indel attribute in the pysam.PileupRead class (0-based query_position)
    MATCH       - 0
    DELETION    - Negative integer on the query position before the deletion for the number of bases
    INSERTION   - Positive integer on the query position before the insertion for the number of bases
    IS_INS      - Returns None-type (no integer)
    IS_DEL      - This program will never iterate over deletions since we are using query positions

    Acceptable cigar strings use MDI - This program does not accepted extended CIGAR format.
    """
    sum = 0
    previous_was_match = False
    for tuple in aligned_segment.cigartuples:
        type = tuple[0]
        length = tuple[1]

        # Valid CIGAR strings alternate between matches and non-matches
        if previous_was_match:
            previous_was_match = False
            if type == 1: # Insertion
                if query_position == (sum - 1):
                    return length
                sum = sum + length
                if query_position < sum:
                    return None
            elif type == 2: # Deletion
                if query_position == (sum - 1):
                    return length * -1
            else:
               raise ValueError("This is not an acceptable CIGAR string: %r" % aligned_segment.cigarstring)
        else:
            if type == 0: # Match
                previous_was_match = True
                sum = sum + length
                if query_position < (sum - 1):
                    return 0
            else:
                raise ValueError("This is not an acceptable CIGAR string: %r, %r" % (aligned_segment.cigarstring, type))
    # Don't forget the last base if it was a match
    if previous_was_match and query_position <= sum:
        return 0
    return IndexError("Query Position %r is out of bounds.")