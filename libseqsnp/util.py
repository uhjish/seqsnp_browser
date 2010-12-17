import sys
import os
import os.path
import shlex
import subprocess


def get_suffix_array_files( path ):
    (directory, prefix) = os.path.split( path )
    suff_list = []
    for f in os.listdir( directory ):
        if f.startswith( prefix ) and f.endswith(".suf"):
            suff_list.append( f.replace(".suf","") )
    return (directory, suff_list)

def get_read_length( directory, suff_list ):
    token = "AAAAAAAAAA"
    idx = directory + "/" + suff_list[0]
    args = ['geoseq', '-p', 'match', '-w', str(len(token)), '-s', token, '-d', idx]    
    pr = subprocess.Popen( args, stdout=subprocess.PIPE )
    #pr.wait()
    seq_lens={}
    for line in pr.stdout:
        if line.startswith(">"):
            continue
        l = len(line.strip())
        if not seq_lens.has_key(l):
            seq_lens[l] = 1
        else:
            seq_lens[l] += 1
    mode = 0
    length = 0
    for k,v in seq_lens.items():
        if v > mode:
            mode = v
            length = k
    return length

