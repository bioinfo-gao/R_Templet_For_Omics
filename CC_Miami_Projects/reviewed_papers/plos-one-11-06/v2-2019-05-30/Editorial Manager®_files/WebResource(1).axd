function getWinWidth(winObj)
{
    if (winObj == null) winObj = window;
    var width = 0;
    if (winObj.innerWidth)
        width = winObj.innerWidth;
    else if (winObj.document.documentElement && winObj.document.documentElement.clientWidth)
        width = winObj.document.documentElement.clientWidth;
    else if (winObj.document.body && winObj.document.body.clientWidth)
        width = winObj.document.body.clientWidth;

    return width;
}
function getWinHeight(winObj)
{
    if (winObj == null) winObj = window;
    var height = 0;
    if (winObj.innerHeight)
        height = winObj.innerHeight;
    else if (winObj.document.documentElement && winObj.document.documentElement.clientHeight)
        height = winObj.document.documentElement.clientHeight;
    else if (winObj.document.body && winObj.document.body.clientHeight)
        height = winObj.document.body.clientHeight;

    return height;
}

function openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct, addWidthHeightToUrl)
{
	// Previously "openCenterWinPct" in clientutil.js
	
	//Popup Arguments
	var scrollArg = (scroll == false) ? '' : ',scrollbars=1';
	var resizeArg = (resize == false) ? '' : ',resizable=1';
	var menubarArg = (menubar == false) ? '' : ',menubar=1';
	var locationArg = (location == false) ? '' : ',location=1';
	var statusArg = (status == false) ? '' : ',status=1';
	var toolbarArg = (toolbar == false) ? '' : ',toolbar=1';
	var minW = minWidth; //Minimum popup width
	var minH = minHeight; //Minimum popup height

	var minW_hit = false; //Is the popup width smaller than the minimum allowed
	var minH_hit = false; //Is the popup height smaller than the minimum allowed
	var winW = 0; //Popup Width
	var winH = 0; //Popup Height
	var winL = 0; //Popup Left Coordinates
	var winT = 0; //Popup Top Coordinates
	var availW = 0; //Width Available to Popup
	var availH = 0; //Height Available to Popup

	//Tabbed browsers (place userAgents here)
	if (navigator.userAgent.indexOf('Opera') != -1)
	{
		var minW = minWidth * 0.9; //Minimum popup width
		var minH = minHeight * 0.8; //Minimum popup height

		//Check if popup is being called from a popup
		if (window.opener == null)
		{
		    availW = getWinWidth();     // Use innerWidth for tabbed browsers
		    availH = getWinHeight();    // Use innerHeight for tabbed browsers
		}
		else
		{
		    availW = getWinWidth(window.opener);    // Use innerWidth for tabbed browsers
		    availH = getWinHeight(window.opener);   // Use innerHeight for tabbed browsers
		}

		//The browsers available browsing dimensions are unknown
		//so if nessassary, adjust minimums to be 70px less than maximum
		if (minW > availW - 70) {minW = availW - 70}
		if (minH > availH - 70) {minH = availH - 70}
	}
	//Regular Browsers
	else
	{
		availW = screen.availWidth; //Use availWidth for regular browsers
		availH = screen.availHeight; //Use availHeight for regular browsers
	}

	//Common Functionality
	winW = availW * pct; //Use Screen Width Percentage
	winH = availH * pct; //Use Screen Height Percentage

	if (winW < minW) {winW = minW; minW_hit = true;}
	if (winH < minH) {winH = minH; minH_hit = true;}

	//Set Window Left
	if (minW_hit) //If popup width will be less than minimum then use the minimum to center
		winL = (availW - minW)/2;
	else
		winL = (availW - winW)/2;

	//Set Window Top
	if (minH_hit) //If popup height will be less than minimum then use the minimum to center
		winT = (availH - minH)/2;
	else
		winT = (availH - winH)/2;

	//Get journal name from the url and add it to the window name. This will make sure there is one and only one
	//window type per journal. This also avoids the re-login problem in the popup Bug 14048 LVL 01242006
	//var the_url = window.location.href;
	//var uri_components = the_url.split("/");
	//var journal_name = uri_components[3];
	//winName = winName + journal_name;

    //Open the popup
	winW = Math.round(winW);
	winH = Math.round(winH);
	if (addWidthHeightToUrl && url != "")
	{
	    url = url + "&width=" + winW + "&height=" + winH;
	}
	var Win = window.open(url,winName,"width="+winW+",height="+winH+",top="+winT+",left="+winL+scrollArg+resizeArg+menubarArg+locationArg+statusArg+toolbarArg);
	if (Win != null) {  // #27696 Alex 20140923 - IE9 if url belongs to a different zone with a different Protected Mode setting value,
	    Win.focus();    // the popup will be opened but the result of the call will be null.
	}
	return Win;
}

function openCenterWinPercentage(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct)
{
    openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct);
}


function openCenterWinStandardRet(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar)
{
	// Previously "openCenterWin2" in clientutil.js
	
	//Popup Arguments
	var scrollArg = (scroll == false) ? '' : ',scrollbars=1';
	var resizeArg = (resize == false) ? '' : ',resizable=1';
	var menubarArg = (menubar == false) ? '' : ',menubar=1';
	var locationArg = (location == false) ? '' : ',location=1';
	var statusArg = (status == false) ? '' : ',status=1';
	var toolbarArg = (toolbar == false) ? '' : ',toolbar=1';

	var winL = 0; //Popup Left Coordinates
	var winT = 0; //Popup Top Coordinates
	var availW = 0; //Width Available to Popup
	var availH = 0; //Height Available to Popup

	//Tabbed browsers (place userAgents here)
	if (navigator.userAgent.indexOf('Opera') != -1)
	{
		//Check if popup is being called from a popup
	    if (window.opener == null)
	    {
	        availW = getWinWidth();     // Use innerWidth for tabbed browsers
	        availH = getWinHeight();    // Use innerHeight for tabbed browsers
	    }
	    else
	    {
	        availW = getWinWidth(window.opener);    // Use innerWidth for tabbed browsers
	        availH = getWinHeight(window.opener);   // Use innerHeight for tabbed browsers
	    }
	}
	//Regular Browsers
	else
	{
		availW = screen.availWidth; //Use availWidth for regular browsers
		availH = screen.availHeight; //Use availHeight for regular browsers
	}

	//Set Window Left
	winL = (availW - winW)/2;

	//Set Window Top
	winT = (availH - winH)/2;

	//Get journal name from the url and add it to the window name. This will make sure there is one and only one
	//window type per journal. This also avoids the re-login problem in the popup Bug 14048 LVL 01242006
	//var the_url = window.location.href;
	//var uri_components = the_url.split("/");
	//var journal_name = uri_components[3];
	//winName = winName + journal_name;

    //Open the popup
	var Win = window.open(url,winName,"width="+winW+",height="+winH+",top="+winT+",left="+winL+scrollArg+resizeArg+menubarArg+locationArg+statusArg+toolbarArg);
	if (Win != null) {  // #27696 Alex 20140923 - IE9 if url belongs to a different zone with a different Protected Mode setting value,
	    Win.focus();    // the popup will be opened but the result of the call will be null.
	}

	// Return window reference
	return Win;
}

function openCenterWinStandard(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar)
{
    openCenterWinStandardRet(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar);
}

function openCenterWin75Percent(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	// Previously "openCenterWin" in clientutil.js
	
	openCenterWinPercentage(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.75);
}

function openCenterWin75PercentRet(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	// Previously "openCenterWin" in clientutil.js
	
	return openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.75);
}

function openCenterWin90Percent(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	// Previously "openCenterWin3" in clientutil.js

	openCenterWinPercentage(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.90);
}

function openCenterWin90PercentRet(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	// Previously "openCenterWin3" in clientutil.js

	return openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.90);
}

function openCenterWin(url, winName, scroll, resize, menubar, location, status, toolbar)
{
    // Previously "openCenterWin" in clientutil.js

    openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.75);
}

function openCenterWinPct(url, winName, scroll, resize, menubar, location, status, toolbar, minWidth, minHeight, pct)
{
    // Previously "openCenterWinPct" in clientutil.js
    openCenterWinPercentageRet(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct);
}



function OpenWindowFixedSize(url, name, width, height, scrollbar) {
    var options = "width=" + width;
    options += ",height=" + height;
    options += ",scrollbars=" + scrollbar;
    options += ",location=no,toolbar=no, menubar=no, resizable=no";

    window.open(url, name, options);
}   
  
function ViewAuthorCorrHistoryPopUp(url)
{
	var aNewWindow = window.open(url, "ViewAuthorCorrHistory","resizable=yes,scrollbars=yes,height=500,width=725");
	
	var width = 725;
	var height = 500;
	
	// Move the popup to the center of the screen
	var x = (screen.availWidth - width) / 2; 
	var y = (screen.avaiHeight - height) / 2 

	aNewWindow.moveTo(x,y);
	
	aNewWindow.focus();
}
