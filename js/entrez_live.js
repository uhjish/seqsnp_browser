
var baseESearch = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
var baseESummary = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
var baseEFetch = "fetchSequence.py";

var MAX_HITS_RETURNED = 100;
var MAX_DROPDOWN_HEIGHT = "100px";
var speciesField;
var queryTextField;
var sequenceField;

function getFieldValue( fld ){
    var res = document.getElementById( fld ).value;
    return(res);
}
function setFieldValue( fld, val ){
    document.getElementById( fld ).value = val;
}

function initializeEntrezLive( spFld, qtFld, sFld ){
    speciesField = spFld;
    queryTextField = qtFld;
    sequenceField = sFld;
} 

function pushToHashSet( hash, id, key, value ){
    if (! id in hash){
        hash[id]=new Array();
    }
    hash[id][key]=value;
}

function loadFasta( chosenId ){
    var url = baseEFetch;
    var params = ({ "gbid":chosenId});
    jQuery.ajax({   type: "GET",
                    url: baseEFetch,
                    async: false,
                    dataType: "text",
                    data: params,
                    success: function( res, st, req ){
                        //window.alert(res);
                        setFieldValue(sequenceField,  res.trim());
                    }, 
                    error: function( req, msg, err){
                        window.alert( msg + err);
                    }
                } );
}

function getMatchingNuclRecords( ){
    var qterm = "("+ getFieldValue( queryTextField ) +"*[Gene Name] OR "+getFieldValue(queryTextField) + "*[Accession])";
    var idList = Array();
    var ssnId = null;
    var qId = null;
    var sp = getFieldValue( speciesField );
    if (sp != "any"){ 
        qterm += " AND " + getFieldValue( speciesField ) + "[orgn]";
    }
   // window.alert(qterm);
    var params = ({ db:"nuccore",
                    term: qterm,
                    retmax: MAX_HITS_RETURNED,
                    usehistory: "y" });
    jQuery.ajax({   type: "GET",
                    url: baseESearch,
                    async: false,
                    dataType: "xml",
                    data: params,
                    success: function( res, st, req ){
                        jQuery(res).find("Id").each( function(){
                                idList.push( jQuery(this).text() );
                        });
                        ssnId = jQuery(res).find("WebEnv").text();
                        qId = jQuery(res).find("query_key").text();
                        //window.alert(ssnId);
                    }, 
                    error: function( req, msg, err){
                        window.alert( msg + err);
                    }
                } );
    //window.alert(ssnId);
    var descs = new Array();
    params = ({ db:"nuccore",
                    WebEnv: ssnId,
                    query_key: "1",
                    retmax: MAX_HITS_RETURNED,
                    rettype: "xml",
                    usehistory: "y" });
    jQuery.ajax({   type: "GET",
                    url: baseESummary,
                    async: false,
                    dataType: "html",
                    data: params,
                    success: function( res, st, req ){
                        var alCt =0;
                        jQuery(res).find("DocSum").each( function(){
                                var id = jQuery(this).find("Id").text();
                                var item = new Array();
                                item["id"]=id;
                                jQuery(this).find("Item").each( function(){
                                    var prop = jQuery(this).attr("Name");
                                    var val = jQuery(this).text();
                                    item[prop]=val;
                                });
                                descs.push(item);
                                
                        });
                    }, 
                    error: function( req, msg, err){
                        window.alert( msg + err);
                    }
                } );

    //return( {"ids": idList, "session":ssnId} );    
    return(descs);                    
}                    

function fetchNuclDetails( srchResults ){

    

}    
    
