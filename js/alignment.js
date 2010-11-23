
var position;

function reverse_complement( s ){
    complement_dna = {A:'T', T:'A', C:'G', G:'C', a:'t', t:'a', c:'g', g:'c', N:'N', n:'n' };
    out = '';
    for (i=s.length; i--;){
        nch = complement_dna[s[i]];
        out += complement_dna[s[i]];
    }
    return out;
}
function formatString(format, str){
    return format.replace("TOKEN",str);
}
function formatReference( pos, refseq){
    base_start=pos-1;
    base_end=pos;
    //format string for snp base matching reference
    fmt_snp = "<font color=green><b>TOKEN</b></font>";
    fmt_flank = "<font color=blue>TOKEN</font>";
    lhighlightMax = 10;
    rhighlightMax = 11;
    lhighlight = position-lhighlightMax;
    if (lhighlight < 0){
        lhighlight =0;
    }
    rhighlight = position+rhighlightMax;
    if (rhighlight > refseq.length){
        rhighlight = refseq.length;
    }
    lBare = refseq.substring(0,lhighlight);
    lFlank = formatString( fmt_flank, refseq.substring(lhighlight,pos) );
    snp = formatString( fmt_snp, refseq.charAt(pos) );
    rFlank = formatString( fmt_flank, refseq.substring(pos+1,rhighlight) );
    rBare = refseq.substring(rhighlight, refseq.length);
    return (lBare + lFlank + snp + rFlank + rBare);
}
function Alignment (qname, start, end, strand, qseq, ref){
        this.strand="+";
        this.qseq = qseq;
        if (strand < 0){
            this.strand = "-";
        }
        this.tstart = start;
        this.tend = end;
        this.qname=qname;
        this.ref=ref;
        this.reflen = ref.length;

    this.formatHighlights = function(pos, qline){
        var base_start=pos;
        var base_end=pos+1;
        var fmt_snp = "<font color=green><b>TOKEN</b></font>";
        if (qline[pos] != this.ref[pos]){
            fmt_snp = "<font color=red><b>TOKEN</b></font>";
        }else{
            var cpos = pos+1;
            while (cpos < qline.length && qline[cpos].toUpperCase() != qline[cpos]){
                base_end += 1;
                cpos+=1;
            }
        }
        var fmt_flank = "<font color=blue>TOKEN</font>";
        var lhighlightMax = 10;
        var rhighlightMax = 11;
        var lhighlight = pos-lhighlightMax;
        if (lhighlight < 0){
            lhighlight =0;
        }
        var rhighlight = pos+rhighlightMax;
        if (rhighlight > this.reflen){
            rhighlight = this.reflen;
        }
        var prefix = "";
        var suffix = "";
        if (this.tstart-1 <= pos && this.tend >= pos){
            prefix = "<a href=# title=\""+this.qname+", "+this.strand+this.tstart+"\" onClick=\"return false;\">";
            suffix = "</a>";
        }
        var lBare = qline.substring(0,lhighlight);
        var lFlank = formatString( fmt_flank, qline.substring(lhighlight,base_start) );
        var snp = formatString( fmt_snp, qline.substring(base_start,base_end) );
        var rFlank = formatString( fmt_flank, qline.substring(base_end,rhighlight) );
        var rBare = qline.substring(rhighlight, qline.length);
        return lBare + prefix + lFlank + snp + rFlank + suffix + rBare;
    }
    this.printHighlighted = function( pos ){
        var qseq = this.qseq;
        var pref = this.tstart-1;
        var suff = this.reflen - this.tend;
        var p_pad = new Array( pref+1 ).join( "-" );
        var s_pad = new Array( suff+1 ).join( "-" );
        var qline = (p_pad+qseq+s_pad).substring(0,this.reflen);
        var retval = this.formatHighlights( pos, qline )
        return retval;
    } 
}

 
function cmpAligns( a, b){
    if (a.tstart == b.tstart){
        return 0;
    } 
    if (a.tstart > b.tstart){
        return 1;
    }else{
        return -1;
    }
}
var MIN_ALIGN_FRACTION=0.9;

function showAlignment( ref, hits, pos ){
    position = pos-1;
    var CLUSTAL_HEADER="<html><head><script language=JavaScript>function do_it(){window.scrollTo("+(position*6.3)+",0);}</script></head><body id=\"aln_body\" onload=do_it()><pre>";
    var CLUSTAL_FOOTER="</pre></body></html>";
    var ref_header = ref["id"].replace(":","_").substring(1,30);
    var refseq = ref["seq"];
    var aligns = new Array();
    for (hidx in hits){
        var h = hits[hidx];
        aligns.push(new Alignment( h["qid"], h["start"], h["end"],h["strand"], h["qseq"], refseq ));
    }
    aligns = aligns.sort(cmpAligns);
    var retStr = CLUSTAL_HEADER;
    retStr += sprintf("\n\n%-20s%-70s%20s\n\n", ref_header.substring(0,15)+"TOP",formatReference(position,refseq), "TOP"+ref_header.substring(0,15));
    for (aidx in aligns){
        var aln = aligns[aidx];
        out =(aln.printHighlighted(position));
        retStr += sprintf("%-20s%-70s%20s\n" ,aln.qname,out, aln.qname);
    }
    retStr += sprintf("\n\n%-20s%-70s%20s", ref_header.substring(0,15)+"BOT",formatReference(position,refseq), "BOT"+ref_header.substring(0,15));
    retStr += CLUSTAL_FOOTER;
    return retStr;
}

