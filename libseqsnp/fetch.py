import urllib
import urllib2
import json
import random

from settings import *

def callSeqSNPService( params ):
    server = random.choice( SEQSNP_SERVERS )
    url = "http://%s/%s" % (server, SEQSNP_ADDRESS)
    data = urllib.urlencode(params)
    req = urllib2.Request( url, data )
    response = urllib2.urlopen(req)
    return json.loads(response.read())
