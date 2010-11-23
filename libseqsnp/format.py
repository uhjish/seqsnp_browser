
def get_ref_flank( refseq, start, end, flank ):
    left = start-flank
    if left < 0:
        left = 0
    right = end + flank
    if right > len(refseq):
        right = len(refseq)
    ret = refseq[left:start] + "[" + refseq[start:end] + "]" + refseq[end:right]
    return ret

def get_snp_rows( ref, snps ):
    snp_table = []
    for s in snps:
        m_frac = "%.2f" % (float(s["mutCount"])/float(s["refCount"]+s["mutCount"]))
        mct_minus = s["mutCount"]-s["mutCountPlusStrand"]
        skew = "%.2f" % (float(s["mutCountPlusStrand"]-mct_minus)/float(s["mutCount"]))
        snpid = "%d:%s>%s" % (s["end"], s["refAllele"], s["mutAllele"].upper() )
        snpid_lnk = "<a href=# onClick=\"zoomToAnnotation(\'%s\');\">%s</a>&nbsp;<font size=-2><a href=# onClick=\"showAlignmentWindow(%d);\">[align]</a></font>" % (snpid, snpid, s["end"])
        snp_table.append( [ snpid_lnk,
                            s["start"],
                            s["end"], 
                            s["refAllele"],
                            s["mutAllele"].upper(),
                            s["refCount"]+s["mutCount"],
                            float(m_frac),
                            float(skew),
                            get_ref_flank( ref["seq"], s["start"], s["end"], 10 ) ] )
    return snp_table
                       
def get_snp_cols():
    names = [ "id", "start", "end", 
            "ref", "mut", "coverage", 
            "mutFreq", "skew", "flank" ]
    cols = []
    for n in names:
        cols.append({"sTitle": n})
    return cols

