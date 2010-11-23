#!/usr/bin/env python
import json

from libseqsnp.fetch import *
from libseqsnp.format import *
from libseqsnp.bp2html import *
from libseqsnp.coverage import *
from mako.template import Template
from mako.lookup import TemplateLookup
from mako import exceptions
import os

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

page = get_required_var("page",form)

mylookup = TemplateLookup(directories=['./tmpl'])


if page == "query":
    print "Content-type: text/html\n\n"
    serve_template('query.tmpl')
elif page == "result":
    print "Content-type: text/html\n\n"
    collection = get_required_var("collection", form)
    basePath = get_required_var("basePath",form)
    seq = get_required_var("seq",form)
    args = {"collection": collection, "basePath": basePath, "seq": seq, "format":"json" }
    calls = callSeqSNPService(args)
    cvg_data = generate_coverage_map( calls[0]["ref"], calls[0]["hits"],calls[0]["filteredMutations"] )
    snp_cols = get_snp_cols()
    snp_rows = get_snp_rows( calls[0]["ref"], calls[0]["filteredMutations"] )
    #print showAlignment( calls[0]["ref"], calls[0]["hits"], 262)
    serve_template('result.tmpl', 
                    cvg_data=json.dumps(cvg_data),
                    seqsnp_data=json.dumps(calls[0]),
                    snp_columns=json.dumps(snp_cols),
                    snp_rows=json.dumps(snp_rows) )         
else:
    raise Exception("Unsupported request for page: %s" % page)
 
