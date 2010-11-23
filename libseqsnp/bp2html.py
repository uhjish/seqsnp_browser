import sys
import math



def reverse_complement( s ):
    complement_dna = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "N":"N", "n":"n" }
    reversed_s = []
    for i in s:
        if i in complement_dna:
            newChar =complement_dna[i]
        else:
            newChar = i
        reversed_s.append(newChar)
    reversed_s.reverse()
    return "".join( reversed_s )

def formatString(format, string):
    return format.replace("TOKEN",string)
def formatReference( pos, refseq):
    base_start=pos-1
    base_end=pos
    #
    #format string for snp base matching reference
    fmt_snp = "<font color=green><b>TOKEN</b></font>"
    fmt_flank = "<font color=blue>TOKEN</font>"
    lhighlightMax = 10
    rhighlightMax = 11
    lhighlight = position-lhighlightMax
    if lhighlight < 0:
        lhighlight =0
    rhighlight = position+rhighlightMax
    if rhighlight > len(refseq):
        rhighlight = len(refseq)
    lBare = refseq[:lhighlight]
    lFlank = formatString( fmt_flank, refseq[lhighlight:pos] )
    snp = formatString( fmt_snp, refseq[pos] )
    rFlank = formatString( fmt_flank, refseq[pos+1:rhighlight] )
    rBare = refseq[rhighlight:]
    return lBare + lFlank + snp + rFlank + rBare

class Alignment:
    def __init__(self, qname, start, end, strand, qseq, ref):
        self.strand="+"
        if strand < 0:
            self.strand = "-"
        qend = None
        self.qseq = qseq
        self.tstart = start
        self.tend = end
        self.qname=qname
        self.ref=ref
        self.reflen = len(ref)
    def formatHighlights( self, pos, qline):
        base_start=pos
        base_end=pos+1
        #
        #format string for snp base matching reference
        fmt_snp = "<font color=green><b>TOKEN</b></font>"
        if qline[pos] != self.ref[pos]:
            #mismatch snp
            fmt_snp = "<font color=red><b>TOKEN</b></font>"
            #check for insertions
        else:
            cpos = pos+1
            while cpos < len(qline) and qline[cpos].upper() != qline[cpos]:
                base_end += 1
                cpos+=1
        #
        fmt_flank = "<font color=blue>TOKEN</font>"
        lhighlightMax = 10
        rhighlightMax = 11
        lhighlight = pos-lhighlightMax
        if lhighlight < 0:
            lhighlight =0
        rhighlight = pos+rhighlightMax
        if rhighlight > self.reflen:
            rhighlight = self.reflen
        prefix = ""
        suffix = ""
        if self.tstart-1 <= pos and self.tend >= pos:
            prefix = "<a href=# title=\"%s, %s\" onClick=\"return false;\">" % ( self.qname, self.strand)
            suffix = "</a>"
        lBare = qline[:lhighlight]
        lFlank = formatString( fmt_flank, qline[lhighlight:base_start] )
        snp = formatString( fmt_snp, qline[base_start:base_end] )
        rFlank = formatString( fmt_flank, qline[base_end:rhighlight] )
        rBare = qline[rhighlight:]
        return lBare + prefix + lFlank + snp + rFlank + suffix + rBare
    def printHighlighted( self, pos ):
        qseq = self.qseq
        #if self.strand =="-":
        #    qseq = reverse_complement(self.qseq)
        pref = self.tstart-1
        suff = self.reflen - self.tend
        #print "pref: %d suff: %d qseq: %s" % (pref,suff,self.qseq)
        qline = (("-"*pref)+qseq+("-"*suff))[:self.reflen]
        retval = self.formatHighlights( pos, qline)
        return retval 
    def __repr__(self):
        return ("seq: %s\n target st: %d en: %d" % (self.qseq, self.tstart, self.tend) )   
        
MIN_ALIGN_FRACTION=0.9
position = None
def showAlignment( ref, hits, pos ):
    global position
    position = pos
    CLUSTAL_HEADER="<html><head><script language=JavaScript>function do_it(){window.scrollTo(%d,0);}</script></head><body onload=do_it()><pre>" % ( position*6.3)
    CLUSTAL_FOOTER="</pre></body></html>"
    ref_header = ref["id"].replace(":","_")[1:30]
    refseq = ref["seq"]
    aligns = {}
    for h in hits:
        #if abs(float(h["end"]-h["start"])) > float(qlen) * MIN_ALIGN_FRACTION:
        aligns[h["qid"]] = Alignment( h["qid"], h["start"], h["end"],h["strand"], h["qseq"], refseq )
    aligns = sorted(aligns.values(), key=lambda aln: aln.tstart)
    retStr = CLUSTAL_HEADER 
    retStr += "\n\n%-20s%-70s\n\n" % (ref_header[0:15]+"TOP",formatReference(position,refseq))
    for align in aligns:
        out =(align.printHighlighted(position)).strip()
        retStr += "%-20s%-70s\n" %(align.qname,out)
    retStr += "\n\n%-20s%-70s\n" % (ref_header[0:15]+"BOT",formatReference(position,refseq))
    retStr += CLUSTAL_FOOTER
    return retStr
