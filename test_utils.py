from libseqsnp.util import *

slist = ['/mnt/geo_sets/JamesIndexer/indexes/CSHL_ENCODE_Stranded_RNASeq/K562_6340_strand1']

length = get_read_length(slist)

print length
