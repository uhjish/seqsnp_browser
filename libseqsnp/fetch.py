import urllib
import urllib2
import json

def callSeqSNPService( params ):
    url = "http://kilauea.mssm.edu:8080/SeqSNPService/seqsnp"
    data = urllib.urlencode(params)
    req = urllib2.Request( url, data )
    response = urllib2.urlopen(req)
    return json.loads(response.read())
