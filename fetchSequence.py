#!/usr/bin/env python
import json

import os
import urllib
import urllib2
import cgi

import cgitb
cgitb.enable()


# CGI header

def get_required_var(var,form):
    if not form.has_key(var):
        raise Exception("Missing required CGI parameter: %s" % var)
    return form.getvalue(var)

def get_optional_var(var,form):
    if not form.has_key(var):
        return False
    return form.getvalue(var)

def serve_template(templatename, **kwargs):
    mytemplate = mylookup.get_template(templatename)
    print mytemplate.render(**kwargs)


form = cgi.FieldStorage()

gbid = get_required_var("gbid",form)

url = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
dat = { "db":"nucleotide",
        "id": gbid,
        "rettype": "fasta" }

dat = urllib.urlencode(dat)
req = urllib2.Request(url, dat)
resp = urllib2.urlopen(req)

print "Content-type: text/plain\n\n"
print resp.read()

