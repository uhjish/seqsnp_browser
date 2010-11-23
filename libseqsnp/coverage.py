import json

def generate_coverage_map( ref, hits, snps ):
    fwd = len(ref["seq"])*[0]    
    rev = len(ref["seq"])*[0] 
    for h in hits:
        if h["strand"] > 0:
            for i in range(h["start"]-1, h["end"]):
                fwd[i]+=1
        else:
            for i in range(h["start"]-1, h["end"]):
                rev[i]+=1
    anns = {}
    for s in snps:
        st = s["start"]
        en = s["end"]
        mut = "%s>%s" % (s["refAllele"],s["mutAllele"].upper())
        name = "%d:%s" % (en, mut)
        if en == st:
            en+=1
        anns[name] = [st, en]       
    return {    "sequence": ref["seq"],
                "forward": fwd, 
                "reverse": rev,      
                "annotations": anns}
            
        
