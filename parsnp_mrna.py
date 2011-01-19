import sys
import shelve
import traceback
from libgeneius import settings
from libgeneius.mysql import GeneiusDb
from libgeneius.simple_geneius import SimpleGeneius
import libgeneius.whereami as whereami
import libgeneius.lookup as lookup
import libgeneius.sequence as sequence
import libgeneius.translate as translate
import libgeneius.codons as codons
from stats import *
from geocounter import *

usage = "python parsnp_mrna.py snpfile build genome_idx snpShelf_list align_shelf patient"

try:
    snpfile = sys.argv[1]
    build = sys.argv[2]
    genome_idx = sys.argv[3]
    shelf_list = sys.argv[4]
    align_shelf = sys.argv[5]
    patient = sys.argv[6]
except:
    print usage
    sys.exit()
mappings = {}

#name    chr     chr_start       chr_end orien   pos     samp    rat     ref_allele      new_allele      Flanking sequence
#ABCC2   chr10   101544294       101544608       +       153     163     1       A       T-1     CCCTTGGGCT[A]CCTATGGCTC


def parseline(line):
    cols = line.strip().split("\t")
    gene = cols[0]
    chr = cols[1]
    ex_st = int(cols[2])
    ex_en = int(cols[3])
    strand = cols[4]
    rel_pos = cols[5]
    hits = int(cols[6])
    mut_freq = float(cols[7])
    ref_allele = cols[8]
    flankseq = cols[10]
    offs = flankseq.find("[")
    flank = flankseq.replace("[","").replace("]","").upper()
    ref_map = geoseq.getcount(flank)
    mut_alleles = {}
    mut_maps = {}
    for mut in cols[9].split(":"):
        if mut[0]=="-":
            allele = "-"
            freq = float(mut[2:])
        else:
            mvals = mut.split("-")
            allele = mvals[0]
            freq = float(mvals[1])
        if allele == "B" or allele == "D":
            continue
        mut_alleles[allele]=freq
        rep_allele = allele
        if rep_allele == "-":
            rep_allele = ""
        mut_flank = flank[:offs] + rep_allele + flank[offs+1:]
        mut_maps[allele] = geoseq.getcount(mut_flank)
    result = {
                "gene": gene,
                "refseq": chr,
                "cds_st": ex_st,
                "cds_en": ex_en,
                "strand": strand,
                "allele": ref_allele,
                "rel_pos": rel_pos,
                "hits": hits,
                "ref_freq": mut_freq,
                "mut_freqs": mut_alleles,
                "flank": flankseq,
                "bare_flank":flank,
                "ref_map": ref_map,
                "mut_maps": mut_maps,
                "flank_offset": offs, 
                "error": None  }
    return result
def get_position(snp):
    snp_position = None
    rel = map(lambda x: int(x), snp["rel_pos"].split(":"))
    if len(rel) ==1:
        snp["snp_start"] = rel[0]-1
        snp["snp_end"] = rel[0]
    else:
        snp["snp_start"] =  rel[0]-1
        snp["snp_end"] = rel[1]


def get_codon_for_snp(snp):
    cdsSt = snp["cds_st"]
    cdsEn = snp["cds_en"]
    pos = snp["snp_start"]
    aa_pos = None
    frame = None
    codon = None
    aa = None
    if pos < cdsSt:
        codon = "5utr"
    elif pos > cdsEn:
        codon = "3utr"
    else:
        codonSt = snp["flank_offset"]
        frame = (pos-cdsSt) % 3
        aa_pos = (pos-cdsSt)/3
        codonSt -= frame
        codon = snp["bare_flank"][codonSt:codonSt+3]
        aa = codons.translate(codon) 
    snp["offset"]=frame
    snp["codon"]=codon
    snp["aa"]=aa
    snp["aa_pos"]=aa_pos

def get_mutation_effects(snp):
    cdsSt = snp["cds_st"]
    cdsEn = snp["cds_en"]
    get_position(snp)
    pos = snp["snp_start"]
    pos_end = snp["snp_end"]
    get_codon_for_snp(snp)
    mutations = []
    for (newbase, al_freq) in snp["mut_freqs"].items():
        ncod = None
        naa = None
        mut_freq = float(al_freq)*float(snp["ref_freq"])
        confs = getSNPpvals(float(snp["hits"]),mut_freq)
        (c_het,p_homo) = confs
        if pos < cdsSt:
            effect = ["5utr"]
        elif pos > cdsEn:
            effect = ["3utr"]
        else:
            if newbase == "-" and (pos_end -pos)%3 == 0:
                effect = ["deletion", "inframe"]
            if newbase == "-" and (pos_end -pos)%3 != 0:
                effect = ["deletion", "frameshift"]
            elif len(newbase) > 1 and (len(newbase)-1) % 3 != 0:
                effect = ["insertion", "frameshift"]
            elif len(newbase) > 1 and (len(newbase)-1) % 3 == 0:
                effect = ["insertion", "inframe"]
            elif newbase == snp["allele"]:
                effect = ["reference"]
                ncod = snp["codon"]
                naa = snpl["aa"]
            elif len(newbase) == 1:
                effect = []
                ofs = snp["offset"]
                ocod = snp["codon"]
                if ocod.startswith("utr"):
                    effect.append("untranslated")
                    ncod=ocod
                    naa = snp["aa"]
                else:
                    ncod = ocod[0:ofs] + newbase + ocod[ofs+1:3]
                    naa = codons.translate(ncod)
                    if naa == snp["aa"]:
                        effect.append("synonymous")
                    elif snp["aa"] == "*":
                        effect.append("runon")
                    elif naa == "*":
                        effect.append("nonsense")
                    else:
                        effect.append("missense")
        mutations.append( { "allele" : newbase,
                            "codon" : ncod,
                            "aa" : naa,
                            "effect" : effect,
                            "conf_het": c_het,
                            "p_homo": p_homo,
                            "frequency": mut_freq,
                            "mut_map": snp["mut_maps"][newbase]  } )
    snp["mutations"] = mutations
     

def get_mappings(snp, aligns):
    maps = aligns[snp["refseq"]]
    coords = []
    for mping in maps:
        coords.append(mping.mrnaToGenomic(snp["snp_start"]))
    return coords



    

geneius_db = GeneiusDb(settings.DB_DATABASE,settings.DB_SERVER,
                       settings.DB_USER,settings.DB_PASSWORD)
genomes_rule = settings.GENOME_PATH

geoseq = GeoTileCounter(genome_idx)

sg = SimpleGeneius()

snp_shelves = []

for line in open(shelf_list, "r"):
    snp_sh = shelve.open(line.strip())
    snp_shelves.append(snp_sh)

aligns = shelve.open(align_shelf)

ct = 0
out = []

for line in open(snpfile,"r"):
    if line.startswith("NO_SNPS_FOUND"):
        break
    if line.startswith("name"):
        print "gene\tregion\tmut_effect\trsid\tref_all\tmut_all\tref_cod\tmut_cod\tref_aa\tmut_aa\ttot_cvg\tmut_freq\tconf_het\tp_homo\tref_map\tmut_map\tflank\trefseq\tmrna_pos\tprot_pos\tgenomic_pos\texon_num\tpileup\tpatient"
        continue
    ct+=1
    snp = parseline(line)
    try:
        get_mutation_effects(snp)
        snp_position = get_mappings(snp,aligns)
    except Exception, pe:
        print >> sys.stderr, "#error: ", line.strip()
        print >> sys.stderr, traceback.format_exc()
        continue 
    for (chrom, strand, gpos,exon) in snp_position:
        for mut in snp["mutations"]:
            snp_start = gpos - 1
            snp_end = snp_start + len(mut["allele"])
            snp_id = "%s:%d-%d" % (chrom,snp_start,snp_end)
            novelty=[]
            for shelf in snp_shelves:
                if shelf.has_key(snp_id):
                    ids = shelf[snp_id]
                    for rsid in ids:
                        if rsid.startswith("rs"):
                            rsid = "<a target=_blank href=http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s>%s</a>" % (rsid,rsid)
                        novelty.append(rsid)
            if not novelty:
                novelty.append("Novel")
            novelty = ",".join(sorted(novelty))
            genomic_link = "<a target=_blank href=http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s:%d-%d>%s:%d-%d:%s</a>" % (chrom,snp_start,snp_end+1,chrom,snp_start,snp_end,strand)
            pileup_link = "<a target=_blank href=../download_sqsh.cgi?ID=%s&refseq=%s&exon=%s&proj=%s&pos=%d>pileup</a>" % (snp["gene"], snp["refseq"], exon, patient ,snp["snp_start"])    
            if snp["codon"] and not snp["codon"].startswith("utr"):
                region = "cds"
            else:
                region = snp["codon"]
            print  "\t".join([
                    snp["gene"],
                    str(region), 
                    "-".join(mut["effect"]),
                    novelty,
                    snp["allele"],
                    mut["allele"],
                    str(snp["codon"]),
                    str(mut["codon"]),
                    str(snp["aa"]),
                    str(mut["aa"]),
                    str(snp["hits"]),
                    ( "%.3f" % float(mut["frequency"])),
                    ( "%.2f" % mut["conf_het"]),
                    ( "%.2e" % mut["p_homo"]),
                    str(snp["ref_map"]),
                    str(mut["mut_map"]),
                    snp["flank"],
                    snp["refseq"],
                    str(snp["snp_start"]),
                    str(snp["aa_pos"]),
                    genomic_link,
                    str(exon),
                    pileup_link,
                    str(patient)
            ])
