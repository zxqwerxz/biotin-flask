# -*- coding: utf-8 -*-
"""Utilities for interacting pysam."""

import os
import pysam

__author__ = 'Jeffrey Zhou'
__copyright__ = 'Copyright (C) 2017, EpigenDx Inc.'
__credits__ = ['Jeffrey Zhou']
__version__ = '0.0.2'
__status__ = 'Development'


###############################################################################
# Wrapper classes around pysam
###############################################################################

class BamFile(pysam.AlignmentFile):
    """Wrapper around the AlignmentFile class.

    Attributes:
        {See the pysam reference for inherited attributes/methods}
        filepath (str): The full path of bam file on disk.

    """

    def __init__(self, filepath, *args, **kwargs):
        """Construct a BAM file."""
        super(BamFile, self).__init__(filename=filepath, *args, **kwargs)
        self.filepath = filepath

    def basename(self):
        """Return the file basename."""
        return os.path.basename(self.filepath)

    def stem(self):
        """Return the file stem, which is the basename w/o extension."""
        basename = os.path.basename(self.filepath)
        return os.path.splitext(basename)[0]

    @classmethod
    def upload(cls, file_h, filename):
        """Upload a SAM/BAM file to disk."""


###############################################################################
# Loading AlignmentFiles
###############################################################################

def load_sams(samfile_list, mode):
    """Load a list samfile paths.

    Parameter:
        samfile_list (list): A list of full paths to sam files.
        mode (str): The mode to open the files.

    Returns:
        A list of AlignmentFile objects.

    """
    result = []
    for samfile in samfile_list:
        result.append(pysam.AlignmentFile(samfile, mode))
    return result


###############################################################################
# Managing Multiple AlignmentFiles
###############################################################################

def get_all_genes(bamlist):
    """Get list of genes from a list of AlignmentFiles.

    Parameter:
        bamlist (list): A list of AlignmentFiles

    Returns:
        A string list of gene names, sorted.

    """
    geneDict = {}
    for bam in bamlist:
        for gene in bam.references:
            if gene not in geneDict:
                geneDict[gene] = True
    return sorted(geneDict.keys())


###############################################################################
# Miscellaneous
###############################################################################

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

def query2refpos(aligned_segment, query_position):
    """
    Returns the reference position at the given query position (0-based).
    Insertions are returned as decimals, for instance:
    152.01 - First insertion after base 152
    152.02 - Second insertion after base 152
    """
    pass

def ref2querypos(aligned_segment, reference_position):
    """
    Returns the reference position at the given query position (0-based).
    My convention for insertions are decimals, for instance:
    152.01 - First insertion after ref base 152
    152.02 - Second insertion after ref base 152
    Return None if there is no base at the given reference position
    """
    # Parse reference positions that are out of range first
    if reference_position < aligned_segment.reference_start or reference_position >= aligned_segment.reference_end:
        return None

    # Reference position is in range
    delta_ref = reference_position - aligned_segment.reference_start
    query_pos = 0
    for tuple in aligned_segment.cigartuples:
        type = tuple[0]
        length = tuple[1]
        if type == 0: # Match
            if delta_ref - (length - 1) <= 0:
                # Need to determine if decimal
                if delta_ref - int(delta_ref) == 0:
                    return query_pos + delta_ref
                else:
                    return None
            delta_ref = delta_ref - (length - 1)
            query_pos = query_pos + length
        elif type == 1: # Insertion
            # Site of interest is only an insertion if delta is a decimal less than 1 but greater than 0
            if delta_ref < 1:
                if round(delta_ref / 0.01, 2) > length:
                    return None
                return int(query_pos - 1 + round(delta_ref / 0.01, 0))
            delta_ref = delta_ref - 1
            query_pos = query_pos + length
        elif type == 2: # Deletion
            if delta_ref - length < 0:
                return None
            delta_ref = delta_ref - length - 1
        else:
            raise ValueError("This is not an acceptable CIGAR string: %r" % aligned_segment.cigarstring)

def get_refbase(aligned_segment, reference_position):
    """
    Returns the base at the given reference position (0-based).
    Return None if the reference_position is out of range.
    Return "-" if there is no base at the given reference_position
    """
    # Parse reference positions that are out of range first
    if reference_position < aligned_segment.reference_start or reference_position >= aligned_segment.reference_end:
        return ""
    pos = ref2querypos(aligned_segment, reference_position)
    if pos is not None:
        return aligned_segment.query_sequence[pos]
    return "-"
