#!/usr/bin/env python
import json

from libseqsnp.fetch import *
from libseqsnp.format import *
from libseqsnp.util import *
from libseqsnp.bp2html import *
from libseqsnp.coverage import *
from mako.template import Template
from mako.lookup import TemplateLookup
from mako import exceptions
import os
import time
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

unbuffered = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stdout = unbuffered

if page == "query":
    print "Content-type: text/html\n\n"
    name = get_optional_var("libname",form)
    path = get_optional_var("libpath", form)
    org = get_optional_var("org", form)
    serve_template('query.tmpl')
    read_length=None
    sys.stdout.flush()
    if path:
        st_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()).strip()
        #find the indices
        (directory, suff_list) = get_suffix_array_files( path )
        suf_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()).strip()
        serve_template('query_loading.tmpl',project=name,suffix_files=suff_list, time_start=st_time, time_rl_start=suf_time) 
        sys.stdout.flush()
        #pull random sequences from the index to get read length
        read_length = getReadLength( path )
        suf_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()).strip()
        serve_template('query_rl.tmpl',read_length=read_length, time_rl_done=suf_time)
        sys.stdout.flush()
    serve_template('query_form.tmpl', libname=name, libpath=path,read_length=read_length)
elif page == "result":
    print "Content-type: text/html\n\n"
    basePath = get_required_var("basePath",form)
    projName = get_optional_var("projName",form)
    if not projName:
        projName = basePath.split("/")[-1]
    seq = get_required_var("seq",form)
    args = {}
    for k in form.keys():
        args[k] = form.getvalue(k)
    res = callSeqSNPService(args)
    try:
        calls = res["results"]
        cvg_data = generate_coverage_map( calls[0]["ref"], calls[0]["hits"],calls[0]["filteredMutations"] )
        snp_cols = get_snp_cols()
        snp_rows = get_snp_rows( calls[0]["ref"], calls[0]["filteredMutations"] )
        #print showAlignment( calls[0]["ref"], calls[0]["hits"], 262)
        serve_template('result.tmpl',
                        projName = projName,
                        cvg_data=json.dumps(cvg_data),
                        seqsnp_data=json.dumps(calls[0]),
                        snp_columns=json.dumps(snp_cols),
                        snp_rows=json.dumps(snp_rows) )        
    except:
        print res 
        raise
else:
    raise Exception("Unsupported request for page: %s" % page)
 
