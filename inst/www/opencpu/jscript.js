//** Global variable **//
//@code collection of cypher
var cypherList = {
    "exactMatch" : "UNWIND keyword AS x WITH x MATCH (node:label) WHERE lower(node.property) = lower(x)",
    "exactCollection" : "UNWIND keyword AS x WITH x MATCH (node:label) WHERE ANY(y IN node.property WHERE lower(y) = lower(x))",
    "regexMatch" : "UNWIND keyword AS x WITH x MATCH (node:label) WHERE lower(node.property) =~ lower(\'.*\'+x+\'.*\')",
    "regexCollection" : "UNWIND keyword AS x WITH x MATCH (node:label) WHERE ANY(y IN node.property WHERE lower(y) =~ lower(\'.*\'+x+\'.*\'))",
    "pathMatch" : ""
    //UNWIND ["coa","co2","Malonyl-CoA"] AS x WITH x MATCH ptw = (n:Metabolite)-[*1..2]-() WHERE lower(n.name) = lower(x) RETURN DISTINCT ptw
    //UNWIND ["coa"] AS x WITH x MATCH ptw = allShortestPaths((n:Metabolite)<-[*]-(nn:Protein)) WHERE lower(n.name) = lower(x) RETURN ptw ORDER BY length(ptw)
}

//@code collection of node properties
var propType = {
    "stringVal" : ["structure","GID","SMILES","equation","IUPAC","description","formula","name","location","InChI","organism","InChIKey","definition"],
    "listVal" : ["type","xref","synonym","describeIn"],
    "booleanVal" : ["manCurate","reversible"],
    "numberVal" : ["km","molWeight"]
}

//@code collection of connection types
var connectList = {
    "biochem" : "",
    "enzcatalyze" : "UNWIND keyword AS x WITH x MATCH (source:Protein {organism:\"species\"}), (target:Metabolite) WHERE lower(target.GID) = lower(x) WITH target, source MATCH ptw = target<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-source RETURN target.GID, target.name, source.GID, source.name, i1 ORDER by source.GID",
    "encgene" : "UNWIND keyword AS x WITH x MATCH (source:Gene {organism:\"species\"}), (target:Metabolite) WHERE lower(target.GID) = lower(x) WITH target, source MATCH ptw = target<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-(i2)<-[:ENCODE]-source RETURN target.GID, target.name, source.GID, source.name, i1, i2 ORDER by source.GID",
    "biopathway" : "UNWIND keyword AS x WITH x MATCH (source:Pathway {organism:\"species\"}), (target:Metabolite) WHERE lower(target.GID) = lower(x) WITH target, source MATCH ptw = target<-[:TRANSFORM|PRODUCE]-(i1)<-[:HAS]-source RETURN target.GID, target.name, source.GID, source.name, i1 ORDER by source.GID"

    //"UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Protein {organism:organism}), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-from RETURN to.GID, to.name, from.GID, from.name, i1 ORDER by from.GID";
    //"UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Gene {organism:organism}), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-(i2)<-[:ENCODE]-from RETURN to.GID, to.name, from.GID, from.name, i1, i2 ORDER by from.GID";
    //"UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Pathway {organism:organism}), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:HAS]-from RETURN to.GID, to.name, from.GID, from.name, i1 ORDER by from.GID";
}

//@code collection of connection colors
var colorList = {
    "biochem" : "#000",
    "enzcatalyze" : "#FF8000",
    "encgene" : "#2E2EFE",
    "biopathway" : "#FF00FF"
}

//** Function **//
//@code collection of cypher
//@param textInput, list of input keywords
//@param searchBy, search by node property
//@param label, node label
//@param matchtype, exact/regex match
//@param organism
//@return cypher string
//@usage
//UNWIND ["coa","co2","Malonyl-CoA"] AS x WITH x MATCH (n:Metabolite) WHERE lower(n.name) = lower(x) RETURN DISTINCT n
//UNWIND ["KEGG:C00010","KEGG:C00083","KEGG:C00011"] AS x WITH x MATCH (n:Metabolite) WHERE ANY(name IN n.xref WHERE lower(name) = lower(x)) RETURN DISTINCT n
//UNWIND ["coa","co2","Malonyl-CoA"] AS x WITH x MATCH (n:Metabolite) WHERE lower(n.name) =~ lower('.*'+x+'.*') RETURN DISTINCT n
//UNWIND ["KEGG:C00010","KEGG:C0008","KEGG:C00011","e"] AS x WITH x MATCH (n:Metabolite) WHERE ANY(name IN n.xref WHERE lower(name) =~ lower('.*'+x+'.*')) RETURN DISTINCT n
//UNWIND ["Fas2"] AS x WITH x MATCH (node:Protein) WHERE lower(node.name) = lower(x) AND node.organism IN ['Arabidopsis thaliana','Homo sapiens'] RETURN DISTINCT node
//UNWIND ["OXSM","KS"] AS x WITH x MATCH (n:Protein) WHERE ANY(name IN n.synonym WHERE lower(name) = lower(x)) AND n.organism IN ["Homo sapiens","Arabidopsis thaliana"] RETURN DISTINCT n
//UNWIND ["FabF"] AS x WITH x MATCH (node:Protein) WHERE lower(node.name) =~ lower('.*'+x+'.*') AND node.organism IN ['Arabidopsis thaliana','Homo sapiens'] RETURN DISTINCT node
//UNWIND ["a"] AS x WITH x MATCH (n:Protein) WHERE ANY(name IN n.synonym WHERE lower(name) =~ lower('.*'+x+'.*')) AND n.organism IN ["Homo sapiens","Arabidopsis thaliana"] RETURN DISTINCT n
function getCypher(textInput,searchBy,label,matchtype,organism) {
    var query = "";
    var isString = $.inArray(searchBy,propType["stringVal"]);
    
    //format textInput
    textInput = decodeURIComponent(textInput);
    textInput = textInput.replace(/\r\n$/, '');//remove last line
    textInput = '[\"'+textInput+'\"]';//as array
    textInput = textInput.replace(/\r\n/g, '\",\"');//as string
    textInput = textInput.replace(/ "/g, '\"'); //replace end whitespace
    
    if (matchtype == 'true' && isString > 0) {
        query = cypherList["exactMatch"];
    }
    else if (matchtype == 'false' && isString > 0) {
        query = cypherList["regexMatch"]; 
    }
    else if (matchtype == 'true' && isString < 0) {
        query = cypherList["exactCollection"]; 
    }
    else{
        query = cypherList["regexCollection"];
    }

    query = query.replace("keyword",textInput);
    query = query.replace("label",label);
    query = query.replace("property",searchBy);
    
    //with organism
    if (organism.length > 0) {
        organism = organism.replace(/,/g, '\",\"');//as string
        query = query + " AND node.organism IN [\"" + organism + "\"]";
    }
    
    query = query + " RETURN DISTINCT node";
    return query;  
}

//@code force exact match
//@param arg, true/false
function forceExactMatch(id,val){ 
    $(id).find('input[name=checkmatch]').prop( 'checked',true );
    $(id).find('input[name=checkmatch]').prop( 'disabled',val );
}

//@code set match type
function setMatchType(id){
    $(id).find('input[name=matchtype]').val( $(id).find('input[name=checkmatch]').is(':checked') ? 'true' : 'false' );      
}

//@code get URL variables
//@param param, name of form input
//@return form values
//@based http://stackoverflow.com/questions/8460265/
function getURLVars(param)
{
    var vars = [], hash;
    var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
    for(var i = 0; i < hashes.length; i++)
    {
        hash = hashes[i].split('=');

        if($.inArray(hash[0], vars)>-1)
        {
            vars[hash[0]]+=","+hash[1];
        }
        else
        {
            vars.push(hash[0]);
            vars[hash[0]] = hash[1];
        }
    }
    var val = vars[param];
    if(val){val = val.replace(/\+/g, ' ');//replace whitespace
    }
    else{val = "";}
    return val;
}

//@code format relationships outputs in result table
//@param value, row, index of tableData
//@return formatted list
function operateRelationship(value, row, index) {
    var lin ="<ul>";   
    for (var x=0; x < value[0].length; x++) {
            lin += '<a href="node.html?GID='+value[1][x]+'"><li>'+value[0][x]+'</li></a>'+value[2][x]+' '+value[3][x].replace(/\[{}\]|\[\"\"\]/g,'');
    }
    lin += "</ul>";
    return lin;
}

//@code format link outputs in result table
//@param value, row, index of tableData
//@return formatted list                        
function operateLink(value, row, index){
    var lin ="<ul>";
        if(typeof value === 'string'){    
            lin += '<a href="https://www.google.com/search?q='+value+'" target="_blank"><li>'+value+'</li></a>';
        }
        else{
            for (var x=0; x < value.length; x++) {
                lin += '<a href="https://www.google.com/search?q='+value[x]+'" target="_blank"><li>'+value[x]+'</li></a>';
            }
        }
    lin += "</ul>";  
    return lin;
}

//@code format name outputs in result table
//@param value, row, index of tableData
//@return formatted list 
function operateName(value, row, index){
    var lin = '<a href="node.html?GID='+value[1]+'">'+value[0]+'</a>';
    return lin;
}

//@code format synonym outputs in result table
//@param value, row, index of tableData
//@return formatted list 
function operateSynonym(value, row, index){
    if (value === undefined) {
        var lin = null;
    }else{  
        var lin ="<ul>";
            if(typeof value === 'string'){    
                lin += value;
            }
            else{ 
                for (var x=0; x < value.length; x++) {
                    lin += '<li>'+value[x]+'</li></a>';
                }
            }
        lin += "</ul>";
    }
    return lin;
}

//@code format edge color
//@param list of relationsip types
//@return edge color 
function formatEdgeColor(edge){
    for (var x=0; x < edge.length; x++) {       
        var rel = edge[x].data.reltype.split(",");        
        var ind = rel.lastIndexOf("NA") + 1; //index of non-NA rel type
        if (ind > 3) {
            var tmprel = rel.splice(3,1);
            ind = tmprel.lastIndexOf("NA") + 1;
        }       
        var ecol = colorList[rel[ind]];        
        edge[x].data["relcolor"] = ecol;       
    }
    return edge;
}

//@code export txt tab files of node attributes and edgelist
//@param list of nodes and edges
function dataToTAB(nodesData, edgesData){
    //for node attributes
    var ntab = '';
    for(key in nodesData[0].data){       
        ntab += key.replace("name","nodeName") + '\t';        
    }
    ntab = ntab.slice(0, -1);
    ntab += '\r\n';
    for (var x = 0; x < nodesData.length; x++){
        for (var key in nodesData[x].data) {
            if (nodesData[x].data.hasOwnProperty(key)) {
              ntab += nodesData[x].data[key] + '\t';
            }
        }
        ntab = ntab.slice(0, -1);
        ntab += '\r\n';
    }

    //for network
    var etab = '';
    for(key in edgesData[0].data){
        etab += key + '\t';
    }
    etab = etab.slice(0, -1);
    etab += '\r\n';
    for (var x = 0; x < edgesData.length; x++){
        for (var key in edgesData[x].data) {
            if (edgesData[x].data.hasOwnProperty(key)) {
              etab += edgesData[x].data[key] + '\t';
            }
        }
        etab = etab.slice(0, -1);
        etab += '\r\n';
    }
    
    //link to txt tab files
    var nuri = 'data:text/plain,' + escape(ntab);
    var euri = 'data:text/plain,' + escape(etab);
    generateFileLink(nuri,"nodeattributes.txt");
    generateFileLink(euri,"edgelist.txt");
}

//@code generate link to txt tab files
//@param data and filename
function generateFileLink(uri,filename){
    var link = document.createElement("a");
    link.href = uri;
    link.style = "visibility:hidden";
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);  
}

//@format textInput and organism for R function
//@param textInput and organism
//@return array contains formatted textInput and organism; e.g. ["["c00025","c00065"]", ""Homo sapiens""]
function formatTextInput(textInput, organism){
    
    //format textInput form
    textInput = decodeURIComponent(textInput);
    textInput = textInput.replace(/\r\n$/, '');//remove last line
    textInput = '[\"'+textInput+'\"]';//as array
    textInput = textInput.replace(/\r\n/g, '\",\"');//as string
    textInput = textInput.replace(/ "/g, '\"'); //replace end whitespace
    organism = '\"'+organism+'\"';
    
    return [textInput, organism];
}



//@base http://jsfiddle.net/motowilliams/7rL2C/
function JSONToTabConvertor(JSONData,file,ShowLabel) {
    //If JSONData is not an object then JSON.parse will parse the JSON string in an Object
    var arrData = typeof JSONData != 'object' ? JSON.parse(JSONData) : JSONData;
    
    var CSV = '';    
        //This condition will generate the Label/Header
    if (ShowLabel) {
        var row = "";
        
        //This loop will extract the label from 1st index of on array
        for (var index in arrData[0]) {
            
            //Now convert each value to string and comma-seprated
            row += index.replace("name","nodeName") + '\t';
        }

        row = row.slice(0, -1);
        
        //append Label row with line break
        CSV += row + '\r\n';
    }
    //1st loop is to extract each row
    for (var i = 0; i < arrData.length; i++) {
        var row = "";
        
        //2nd loop will extract each column and convert it in string comma-seprated
        for (var index in arrData[i]) {
            row += arrData[i][index] + '\t';
        }

        row = row.slice(0, -1);
        //add a line break after each row
        CSV += row + '\r\n';
    }

    if (CSV == '') {        
        alert("Invalid data");
        return;
    }   
    
    //Initialize file format you want csv or xls
    var uri = 'data:text/csv;charset=utf-8,' + escape(CSV);
    generateFileLink(uri,file);
    
}