def getCigarLen(cigar):
    temp = ""
    seqPos = 0
    for c in cigar:
        if c.isdigit():
            temp += c
        elif c == "M":
            seqPos += int(temp)
            temp = ""
        elif c == "I":
            seqPos += int(temp)
            temp = ""
        elif c == "D":
            temp = ""
    return seqPos

def getBase(read, pos):
    start = read.reference_start + 1
    end = read.reference_end + 1
    cigar = read.cigarstring
    seq = read.query_sequence
    if (pos - start < 0) or (end - pos < 0):
        return ""
    target = pos - start
    temp = ""
    seqPos = 0
    refPos = 0
    for c in cigar:
        if c.isdigit():
            temp += c
        elif c == "M":
            if refPos + int(temp) > target:
                if seqPos+target-refPos < len(seq):
                    return seq[seqPos+target-refPos]
                return ""
            refPos += int(temp)
            seqPos += int(temp)
            temp = ""
        elif c == "I":
            seqPos += int(temp)
            temp = ""
        elif c == "D":
            if refPos + int(temp) > target:
                return "-"
            refPos += int(temp)
            temp = ""
    return ""