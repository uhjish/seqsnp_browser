import urllib
import urllib2
import json
import sys
import random

from settings import *

def callSeqSNPService( params ):
    server = random.choice( SEQSNP_SERVERS )
    url = "http://%s/%s" % (server, SEQSNP_ADDRESS)
    data = urllib.urlencode(params)
    req = urllib2.Request( url, data )
    response = urllib2.urlopen(req)
    resp = response.read()
    try:
        ret = json.loads(resp)
    except:
        raise Exception("Invalid json returned:\n"+ resp)
    return ret

def getReadLength( basePath ):
    server = random.choice( SEQSNP_SERVERS )
    url = "http://%s/%s" % (server, SEQSNP_ADDRESS)
    data = urllib.urlencode({"app":"readlength","basePath":basePath})
    req = urllib2.Request( url, data )
    response = urllib2.urlopen(req)
    resp = response.read()
    try:
        ret = json.loads(resp)
    except:
        raise Exception("Invalid json returned:\n"+ resp)
    return ret

