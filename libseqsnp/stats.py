import scipy.stats
import math
import sys

def getSNPpvals( cvg, freq):
    confidence = scipy.stats.binom.sf(round(cvg*freq), cvg, 0.5)
    if confidence > 0.5:
        confidence = 1-confidence
    conf_het = confidence * 200
    p_homo = confidence
    return (conf_het, p_homo)

def __main__():
    cvg = sys.argv[1]
    freq = sys.argvp[2]
    (v1,v2)=getSNPpvals(cvg,freq)
    print cvg, freq, v1, v2
    
