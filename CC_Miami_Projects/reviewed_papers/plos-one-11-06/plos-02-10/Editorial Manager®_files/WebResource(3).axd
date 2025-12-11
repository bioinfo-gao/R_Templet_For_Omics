

var xmlHttpFactories = [
 function() { return new XMLHttpRequest(); }, // IE 7.0 and everything non-Microsoft 
 function()  { return new ActiveXObject("Msxml3.XMLHTTP");},
 function() { return new ActiveXObject("Msxml2.XMLHTTP.6.0"); },
 function()  { return new ActiveXObject("Msxml2.XMLHTTP.3.0"); },
 function()  { return new ActiveXObject("Msxml2.XMLHTTP"); },
 function()  { return new ActiveXObject("Microsoft.XMLHTTP"); }];
 
 function createXmlHttpObject()
 {
    var xmlHttp = false;
    
    for (var i=0; i < xmlHttpFactories.length; i++)
    {
        try
        {        
             xmlHttp = xmlHttpFactories[i]();
        }
        catch(e)
        {
            continue;
        }
        break;
    }
    
    return xmlHttp;
 
 }
 
 function configureXmlHttpRequest(xmlHttpReqObject, url, postData, asynch)
 {
    var sendMethod = (postData) ? 'POST' : 'GET';
   
    xmlHttpReqObject.open(sendMethod, url, asynch);
    xmlHttpReqObject.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');
     //Spec 14.0-14 Added length check as part of spec to prevent long standing error from occuring . 
    if (postData && postData.length!=null)
    { 
        xmlHttpReqObject.setRequestHeader("Content-length", postData.length);
    }
 //Spec 14.0-14 Removed as part of spec to prevent long standing error from occuring .     xmlHttpReqObject.setRequestHeader("Connection", "close");
    
 }
 
function sendXmlHttpRequest(callbackMethod, url, postData, asynch) 
{
   
    var xmlHttpReqObject = false;
    xmlHttpReqObject = createXmlHttpObject();
            
    if (xmlHttpReqObject)
    {
        configureXmlHttpRequest(xmlHttpReqObject, url, postData, asynch);
        
        xmlHttpReqObject.onreadystatechange = 
		    function() 
		    {
                if (callbackMethod && xmlHttpReqObject.readyState == 4 && xmlHttpReqObject.status == 200 ) 
		        {    
        		    callbackMethod(xmlHttpReqObject);
		        } 
		    };
		    
		xmlHttpReqObject.send(postData);
        
        delete xmlHttpReqObject;
    }
}


function xmlHTTPPost(updateFunction, url, querystring) 
{
    var xmlHttpReqObject = false;
    
    xmlHttpReqObject = createXmlHttpObject();
    if (xmlHttpReqObject)
    {

        configureXmlHttpRequest(xmlHttpReqObject,url,true,true);
        xmlHttpReqObject.onreadystatechange = 
		    function() 
		    {
			    if (xmlHttpReqObject.readyState == 4 && xmlHttpReqObject.status == 200 ) 
			    { 
				    updateFunction(xmlHttpReqObject.responseText, "hello");
			    }
		    }
    	
        xmlHttpReqObject.send(querystring);
                
        delete xmlHttpReqObject;
        
    }
}

function xmlHTTPPostWithContext(updateFunction, url, querystring, context) 
{
    var xmlHttpReqObject = false;
        	
	xmlHttpReqObject = createXmlHttpObject();
    
    if (xmlHttpReqObject)
    {	
        configureXmlHttpRequest(xmlHttpReqObject, url, true,true);
        
        xmlHttpReqObject.onreadystatechange = 
		    function() 
		    {
			    if (xmlHttpReqObject.readyState == 4 && xmlHttpReqObject.status == 200 ) 
			    { 
				    updateFunction(xmlHttpReqObject.responseText, context);
			    }
		    }
    		
        xmlHttpReqObject.send(querystring);
        delete xmlHttpReqObject;
    }
    return true ;
}
