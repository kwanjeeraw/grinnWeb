function buildNetwork(textInput){
    //format textInput|result from convert ID function
    textInput = decodeURIComponent(textInput);
    textInput = textInput.replace(/\r\n$/, '');//remove last line
    textInput = '[\"'+textInput+'\"]';//as array
    textInput = textInput.replace(/\r\n/g, '\",\"');//as string
    textInput = textInput.replace(/ "/g, '\"'); //replace end whitespace
    
    var qprotein = "UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Protein), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-from RETURN to.GID, to.name, from.GID, from.name, i1 ORDER by from.GID";
    var qgene = "UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Gene), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:CATALYZE]-(i2)<-[:ENCODE]-from RETURN to.GID, to.name, from.GID, from.name, i1, i2 ORDER by from.GID";
    var qpathway = "UNWIND [\"met7\",\"met1\"] AS x WITH x MATCH (from:Pathway), (to:Metabolite) WHERE lower(to.GID) = lower(x) WITH to, from MATCH ptw = to<-[:TRANSFORM|PRODUCE]-(i1)<-[:HAS]-from RETURN to.GID, to.name, from.GID, from.name, i1 ORDER by from.GID";

    var req = ocpu.rpc("connectNodes", {querystring: query}, function(object) {});
}