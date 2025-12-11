
var FALSE = 0;
var False = 0;
var TRUE  = 1;
var True  = 1;

var OPTIONAL    = 0;
var REQUIRED    = 1;
var HIDDEN      = 2;

var REVIEWER_SELECTION_MODE_SEARCH = 0;
var REVIEWER_SELECTION_MODE_BY_CLASS = 1;
var REVIEWER_SELECTION_MODE_BY_PERSONAL_CLASS = 2;
var REVIEWER_SELECTION_MODE_SUGGEST = 3;
var REVIEWER_SELECTION_MODE_PREVIOUS = 4;
var REVIEWER_SEARCH_ENTIRE_DATABASE = 999999999;


var REV_INIT_TYPE_ByExistingReviewerCollections = 1;
var REV_INIT_TYPE_ByReviewerRoleIDLists = 2;
var REV_INIT_TYPE_ByXmlString = 3;
var REV_INIT_TYPE_ByPeopleIDList = 4;
var REV_INIT_TYPE_Empty = 5;

var REV_COLLETION_TYPE_InvitedReviewers = 1;
var REV_COLLETION_TYPE_AlternateReviewers = 2;
var REV_COLLETION_TYPE_AlternateLinkedReviewers = 3;
var REV_COLLETION_TYPE_ProposedReviewers = 4;
var REV_COLLETION_TYPE_All = 5;
var REV_COLLETION_TYPE_AssignedNotInvitedReviewers = 6; //9.0-53
		
var REV_TYPE_InvitedReviewer = 1;		
var REV_TYPE_AlternateReviewer = 2;		
var REV_TYPE_AlternateLinkedReviewer = 3;		
var REV_TYPE_ProposedReviewer = 4;		
var REV_TYPE_None = 5;		
var REV_TYPE_AssignedNotInvitedReviewer = 6;	//9.0-53

var REV_PROCESSACTIONTYPE_ChangeSearch = 1;
var REV_PROCESSACTIONTYPE_ConfirmAndSelect = 2;
var REV_PROCESSACTIONTYPE_CancelSearchResults = 3;
var REV_PROCESSACTIONTYPE_AddReviewerRoles = 4;
var REV_PROCESSACTIONTYPE_ChangeSelections = 5;
var REV_PROCESSACTIONTYPE_CancelAssignReviewers = 6;

var REV_PROMOTE_TYPE_PromoteToAlternate = 1;
var REV_PROMOTE_TYPE_PromoteToInvite = 2;
var REV_PROMOTE_TYPE_Unspecified = 3;

var REV_PREF_ShowMatchesOnly = 1;
var REV_PREF_Supress = 2;
var REV_PREF_ShowAll = 3;

var OPERATION_ADD     = 1;
var OPERATION_DELETE  = 2;
var OPERATION_EDIT    = 3;
var OPERATION_REORDER = 4;
var OPERATION_SAVE = 5;
var OPERATION_CLICKED_UP_ARROW   = 1;
var OPERATION_CLICKED_DOWN_ARROW = 2;

var ALTERNATE_OPERATION_REORDER = 1;
var ALTERNATE_OPERATION_ADD     = 2;
var ALTERNATE_OPERATION_REMOVE  = 3;

var UPDATE_INSTRUCTIONS			   = 1;
var REVERT_TO_DEFAULT_INSTRUCTIONS = 2;


var INVITATION_MODE = 0;
var AGREED_MODE = 1;


var AUTHOR_SEARCH_MODE_BY_NAME_ADDR   = 1;
var AUTHOR_SEARCH_MODE_BY_PERS_CLASS  = 2;

var AUTHOR_SEARCH_MODE_BY_NAME_ADDR_ALTERNATE = 3;
var AUTHOR_SEARCH_MODE_BY_PERS_CLASS_ALTERNATE = 4;


var ARTICLE_TYPE_FAMILY_REGULAR  = 1;
var ARTICLE_TYPE_FAMILY_PROPOSAL = 2;
var ARTICLE_TYPE_FAMILY_CONSUB   = 3;	
var STEP_ONE   = 1;
var STEP_TWO   = 2;
var STEP_THREE = 3;
var NOTHING_IS_CHOSEN_VALUE = -1

var RF_AUTHOR = "AUTHOR";
var RF_REVIEWER = "REVIEWER";
var RF_EDITOR = "EDITOR";
var RF_PUBLISHER = "PUBLISHER";

var RF_ID_AUTHOR = 1;
var RF_ID_REVIEWER = 2;
var RF_ID_EDITOR  = 3;
var RF_ID_PUBLISHER = 4;


var CORRES_TYPE_EXTERNAL_EDITORIAL = 1;
var CORRES_TYPE_EXTERNAL_PRODUCTION = 2;
var CORRES_TYPE_INTERNAL_EDITORIAL = 3;
var CORRES_TYPE_INTERNAL_PRODUCTION =4;

var EXT_CORRES_ASSOC_TYPE_SUBMISSION = 1;
var EXT_CORRES_ASSOC_TYPE_SCHEDULEGROUP = 2;

	
	


var PRIOR_TO_SUBMISSION = 0;
var AFTER_AUTHOR_ACCEPTS = 1;

var MAX_DIGITS = 18;
var MIN_DIGITS = 2;



// Custom Image types
var CUSTOM_IMAGE_LETTER = 1;
var CUSTOM_IMAGE_LIVE_LOGIN_THUMBNAIL = 2;
var CUSTOM_IMAGE_WORKING_LOGIN_THUMBNAIL = 3;
var CUSTOM_IMAGE_PDF_WATERMARK = 4; // 13.0-07  20150803 TimH

// PDF Image State/Status
// 13.0-07  20150803 TimH
var PDF_IMAGE_STATE_NONE = 0;
var PDF_IMAGE_STATE_NEW = 1;
var PDF_IMAGE_STATE_COMMITTED = 2;
var PDF_IMAGE_STATE_REMOVED = 3;
var PDF_IMAGE_STATE_SELECTEDFORPDF = 4; 

// Image Orientation
// 13.0-07  20150803 TimH
var PDF_IMAGE_ORIENTATION_PORTRAIT = 1;
var PDF_IMAGE_ORIENTATION_LANDSCAPE = 2;

// PDF Page Type the Watermark Image will
// be placed on.
// 13.0-07  20150803 TimH
var PDF_PAGETYPE_FOR_IMAGE_NONE	= 0;
var PDF_PAGETYPE_FOR_IMAGE_AUTHOREDITOR_COVER = 1;
var PDF_PAGETYPE_FOR_IMAGE_AUTHOREDITOR_PORTRAIT = 2;
var PDF_PAGETYPE_FOR_IMAGE_AUTHOREDITOR_LANDSCAPE = 3;
var PDF_PAGETYPE_FOR_IMAGE_REVIEWER_COVER	= 4;
var PDF_PAGETYPE_FOR_IMAGE_REVIEWER_PORTRAIT = 5;
var PDF_PAGETYPE_FOR_IMAGE_REVIEWER_LANDSCAPE = 6;

// UPLOAD schemes
var UPLOAD_SCHEME_SUBMISSION_FILES = 1;
var UPLOAD_SCHEME_REVIEWER_ATTACHMENT = 2;
var UPLOAD_SCHEME_UPLOAD_WF_COVER_PAGE = 3;
var UPLOAD_SCHEME_SUBSTITUTE_ITEMS = 4;
var UPLOAD_SCHEME_CUSTOM_IMAGE = 5;
var UPLOAD_SCHEME_COMPANION_FILE = 8;		// Spec 4.0-04. 20051103 LEH

// Upload actions
var CUSTOM_IMAGE_UPLOAD = 1;
var CUSTOM_IMAGE_EDIT = 2;
var CUSTOM_IMAGE_REPLACE = 3;

// DOWNLOAD schemes
var DOWNLOAD_SCHEME_DEFAULT = 1;
var DOWNLOAD_SCHEME_SCRATCH_FILES = 2;
var DOWNLOAD_SUBSTITUTE_ITEM = 3;
var DOWNLOAD_SCHEME_CUSTOM_IMAGE = 4;
var DOWNLOAD_SCHEME_WEBFIRST_SUBSCRIBER = 5;
var DOWNLOAD_SCHEME_WEBFIRST_PPV = 6;
var DOWNLOAD_SCHEME_LOGIN_FILE = 7;
var DOWNLOAD_SCHEME_REVIEWER_ATTACHMENT = 8;
var DOWNLOAD_SCHEME_AMENDED_FILES = 9;
var DOWNLOAD_SCHEME_COMPANION_FILE = 10;	// Spec 4.0-04. 20051109 LEH
var DOWNLOAD_SCHEME_ASSIGNMENT_FILE = 11;	// Spec 6.0-23-17 (used for tracking last download)
var DOWNLOAD_SCHEME_SG_ASSIGNMENT_FILE = 12;	// Spec 7.0-32 Used for tracking download for schedule group task assignment files
var DOWNLOAD_SCHEME_LETTERATTACHMENT = 13;      // spec 8.0-20 new for Letter Attachments



function isPosNonZeroInt(inputString, maxNumDigits)
{
    return ((isInteger(inputString, maxNumDigits)) && (inputString > 0));
}

function isPosInt(inputString, maxNumDigits)
{
    return ((isInteger(inputString, maxNumDigits)) && (parseInt(inputString) >= 0));
}

function isInteger(inputString, maxNumDigits)
{
    if (typeof(maxNumDigits) == "undefined")
    {
        
        var regExp = new RegExp("^(\\s*-?\\d+\\s*)$");
    }
    else
    {
        
        var regExp = new RegExp("^\\s*-?\\d\{1," + maxNumDigits + "\}\\s*$");
    }

    return regExp.test(inputString);
}

function isFloat(inputString)
{

    var regExp = new RegExp("^((-?[0-9]+(\\.[0-9]*)?)|(-?([0-9]*\\.)?[0-9]+))$");

    return regExp.test(inputString);
}

function isDate(inputString)
{
    
    var regExp28 = new RegExp("^(\\s*(0?[1-9]|1[012])[\- \/\.](0?[1-9]|1[0-9]|2[0-8])[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExp30 = new RegExp("^(\\s*(0?[13-9]|1[012])[\- \/\.](29|30)[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExp31 = new RegExp("^(\\s*(0?[13578]|1[02])[\- \/\.](31)[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExpLeap = new RegExp("^(\\s*(0?2)[\- \/\.](29)[\- \/\.]([1-9][0-9])?(([0248][048])|([13579][26]))\\s*)$");

    if (regExp28.test(inputString) || regExp30.test(inputString) || regExp31.test(inputString)) 
    {
        return true;
    }
    else if (regExpLeap.test(inputString))
    {
        var mdy = inputString.split(/[\- \/\.]/);
        var month = null;
        var date = null;
        var year = null;
        
        if (mdy && mdy.length == 3)
        {
            month = parseInt(mdy[0].replace(/^\s*0+/,""));
            date = parseInt(mdy[1].replace(/^\s*0+/,""));
            if (month != 2 || date != 29)
            {
                return false;
            }
            
            var year = parseInt(mdy[2].replace(/^\s*0+/,""));
            var today = new Date();
            var currYear = today.getFullYear();
            var currMillenium = Math.floor(currYear / 1000);
            if (year < currMillenium * 1000)
            {
                year = currMillenium * 1000 + (year % 100);
            }
            
            
            if (0 == year % 100)
            {
                return ( 0 == year % 400);
            }
            else
            {
                return true;
            }            
        }
    }
	return false;
}


function isValidStandardDate(dateString, mustBeFuture)
{
	var dateString = "" + dateString;
	dateString = dateString.replace(/\s+/g, "");
	
	if (dateString == "")
	{
		return true;
	}
	
	var patternMatch = dateString.match(/\d{2}\/\d{2}\/\d{4}/);
	
	if (patternMatch == null)
	{
		return false;
	}
	
	var inputDate = new Date(dateString);
			    				    		    						    		    		
	if (inputDate == null)
	{
		return false;
	}
	
	var year = inputDate.getFullYear();
	var month = inputDate.getMonth() + 1;
	var day = inputDate.getDate();
	
	if (isNaN(year) || isNaN(month)|| isNaN(day))
	{
		return false;
	}
	
	var yyyy = '' + year;
	var mm = (month > 9) ? month  : '0' + month;
	var dd = (day > 9) ? day : '0' + day;
	
	var mdy = mm + '/'  + dd + '/' + yyyy;
	
	if (dateString != mdy)
	{
		return false;
    }  				   
	
	if (mustBeFuture)
	{
		var currDate = new Date();
		currDate.setHours(0);
		currDate.setMinutes(0);
		currDate.setSeconds(0);
		currDate.setMilliseconds(0);
		
		return (inputDate >= currDate);
	}
	
	return true;
}


function tryParseDate(inputString)
{
    
    var regExp28 = new RegExp("^(\\s*(0?[1-9]|1[012])[\- \/\.](0?[1-9]|1[0-9]|2[0-8])[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExp30 = new RegExp("^(\\s*(0?[13-9]|1[012])[\- \/\.](29|30)[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExp31 = new RegExp("^(\\s*(0?[13578]|1[02])[\- \/\.](31)[\- \/\.]([1-9][0-9])?[0-9]{2}\\s*)$");
    
    var regExpLeap = new RegExp("^(\\s*(0?2)[\- \/\.](29)[\- \/\.]([1-9][0-9])?(([0248][048])|([13579][26]))\\s*)$");

    var isValidDateFormat = false;
    var isLeapYear = false;
    var dateObject = null;
    
    if (regExp28.test(inputString) || regExp30.test(inputString) || regExp31.test(inputString)) 
    {
        isValidDateFormat = true;
    }

    else if (regExpLeap.test(inputString))
    {
        isValidDateFormat = true;
        isLeapYear = true;
    }
    
    if (isValidDateFormat)
    {
    
        var mdy = inputString.split(/[\- \/\.]/);
        var month = null;
        var date = null;
        var year = null;
            
        if (mdy && mdy.length == 3)
        {
            month = parseInt(mdy[0].replace(/^\s*0+/,""));
            date = parseInt(mdy[1].replace(/^\s*0+/,""));
            
            var year = parseInt(mdy[2].replace(/^\s*0+/,""));
            var today = new Date();
            var currYear = today.getFullYear();
            var currMillenium = Math.floor(currYear / 1000);
            if (year < 100)
            {
                year = currMillenium * 1000 + (year % 100);
            }
                    
            if (!isLeapYear || (0 != year % 100) || ( 0 == year % 400))
            {
                dateObject = new Date(year, month-1, date);
                
                if (isNaN(dateObject))
                    dateObject = null;
            }
              
        }
    }
    
    return dateObject;              
}


function trim(stringToTrim)
{
    if (stringToTrim == null)
        return null;

	stringToTrim = stringToTrim.replace(/^\s+/, '');	
	stringToTrim = stringToTrim.replace(/\s+$/, '');	

    return stringToTrim;
}

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

	function popupPublishTargetInfo(docID, ms_num, sPage, ePage, numPages, tocPos)
	{
		var name = "med_publish_info.asp?docID=" + docID + "&ms_num=" + ms_num + "&sPage=" + sPage +
				   "&ePage=" + ePage + "&numPages=" + numPages + "&tocPos=" + tocPos;
		openCenterWin(name,"publish_information",1,1,0,0,0,0);
	}


function popupZipDownload(QS)
{
    var name = "publisherZip.asp?" + QS;
	openCenterWin(name,"ZipDownload",1,1,0,0,0,0)
}

  

function popupDetailsWindow(docID, ms_num, sectionID)
{
    var name = "EMDetails.aspx?docid=" + docID + "&ms_num=" + ms_num + "&sectionID=" + sectionID;
    openCenterWin(name,"popupDetailsWindow",1,1,0,0,0,0)
}

function popupDetailsWindow2(docID, ms_num, sectionID)
{
    var name = "EMDetails.aspx?docid=" + docID + "&ms_num=" + ms_num + "&sectionID=" + sectionID;
    var pathArray = window.location.pathname.split('/');
    if (pathArray.length > 3)
    {
        for(var i=3; i < pathArray.length; i++)
        {
            name = '../' + name;
        }
    }
    openCenterWin(name,"popupDetailsWindow2",1,1,0,0,0,0)
}

function popupReprintForm(docID)
{
    var name = 'chooseReprintType.asp?docid=' + docID;
	openCenterWin(name,"popupReprintForm",1,1,0,0,0,0)
}

function popupReprintFormFromEM(docID)
{
    var name = 'PM/chooseReprintType.asp?docid=' + docID;
    openCenterWin(name,"popupReprintFormFromEM",1,1,0,0,0,0)
}

function popupReprintHelp(authorReprint)
{
	var name = "";

	if (authorReprint)
		name = "http://www.docurights.com/APrintHelp.html";
	else
		name = "http://www.docurights.com/CPrintHelp.html";

	openCenterWin(name,"popupReprintFormFromEM",1,1,0,0,0,0)
}


function popupReprintPricing(type)
{
    var name = 'showReprintPrices.asp?type=' + type;
    openCenterWin(name,"popupReprintPricing",1,1,0,0,0,0);
}

function popupReviewerInfo(peopleID, docID, jrnlID)
{
    var url = "reviewer_info.asp?rv=R";
    url += "&id=" + peopleID;
    url += "&detail=1";
    url += "&docid=" + docID;

	if (typeof(jrnlID) != "undefined")
		url += "&jrnlID=" + jrnlID;

    openCenterWin2(url,"popupReviewerInfo",750, 500, 1,1,0,0,0,0)
}


function popupReviewDetails(unBlindReviewer, peopleID, docID, revision, eligibleForReOpen)
{
    

    
    var url = "reviewer_comments.asp?docid=" + docID + "&peopleid=" + peopleID +
						"&rev=" + revision + "&unBlindReviewer=" + unBlindReviewer +
                        "&eligibleForReOpen=" + ((eligibleForReOpen != null && eligibleForReOpen != "undefined") ? eligibleForReOpen : "false");

    openCenterWin(url, "reviewer_comments2",1,1,0,0,0,0)

}

function popupRequestUnregisteredReviewer(docID, revision, ms_num)
{
    var url;
    url = "requestUnregReviewer.asp?docid=" + docID;
    url += "&rev=" + revision;
    url += "&ms_num=" + ms_num;

    openCenterWin(url,"popupRequestUnregisteredReviewer",1,1,0,0,0,0);
}

function popupDeclineReason(roleRevuvID)
{
    var name = "reviewerDeclineReason.asp?roleid=" + roleRevuvID;
    openCenterWinPct(name,"reviewerDeclineReason",1,1,0,0,0,0, 400, 200, 0.25);
}

function popupUnavailableReason(peopleID, jrnlID)
{
    var url = "UnavailabilityPopup.aspx?peopleID=" + peopleID;

	if (typeof(jrnlID) != "undefined")
		url += "&jrnlID=" + jrnlID;

    openCenterWinPct(url,"unavailablePopup",1,1,0,0,0,0, 900, 400, 0.25);
}

function popupClassifications(docID, ms_num, sessionThreadParam)
{
    popupClassificationsRet(docID, ms_num, sessionThreadParam); 
}

function popupClassificationsRet(docID, ms_num, sessionThreadParam)
{
    var name = "SelectSubmissionClassifications.aspx";
	if (docID > 0)
	{
		name = name + "?docid=" + docID + "&ms_num=" + ms_num;
	}
	else
	{
		name = name + "?returnSelectedClassif=true";
	}

    name += (sessionThreadParam) ? '&' + sessionThreadParam: "";

	return openCenterWinPctRet(name,"editor_classify",1, 1, 0, 0, 0, 0, 750, 650, 0.85, 1); 
}

function popupHistoryWindow(docID, ms_num, stat, corr)
{

	if (stat == true)
		stat = 1;
	else
		stat = 0;
	if (corr == true)
		corr = 1;
	else
		corr = 0;

    var name = "doc_history.asp?docid=" + docID + "&ms_num=" + ms_num + "&stat=" + stat +
               "&corr=" + corr;
     openCenterWin(name,"doc_history",1,1,0,0,0,0);
}

function popupScheduleGroupHistoryWindow(scheduleGroupID)
{
    var name = "ScheduleGroupHistory.aspx?scheduleGroupID=" + scheduleGroupID;
     openCenterWinPct(name,"scheduleGroupHistory",1,1,0,0,0,0, 750, 500, 0.90);
}

function popupScheduleGroupFileInventory(scheduleGroupID)
{
    var name = "ScheduleGroupFileInventory.aspx?scheduleGroupID=" + scheduleGroupID;
	openCenterWin(name,"scheduleGroupFileInventory",1,1,0,0,0,0);
}

function popupViewReviewerComments(docID, ms_num, revision)
{
    var name = "listReviewersWithReview.asp?docid=" + docID + "&ms_num=" + ms_num + "&rev=" +
               revision;
    openCenterWin(name,"document_details",1,1,0,0,0,0);
}

function viewReviewsAndComments(docID, revision, manuscriptNumber)
{
    var name = "viewReviewerAndEditorComments.asp?docid=" + docID + "&rev=" + revision + "&ms_num=" + encodeURIComponent(manuscriptNumber);
    openCenterWin(name,"viewReviewsAndComments",1,1,0,0,0,0)
}

function viewEditorComments(editorID, docID, revision)
{
    var name = "editor_comments.asp?peopleid=" + editorID + "&docid=" + docID + "&rev=" + revision;
	openCenterWin(name,"editor_comments",1,1,0,0,0,0);
}

function popupAQCFileResults(AQCDataID)
{
    var name = "viewAQCResults.aspx?aqcdataid=" + AQCDataID;
    openCenterWin(name,"viewFileAQC",1,1,0,0,0,0);
}

function popupAQCMSResults(docID, revision, approvalType, refreshPage)
{
    var name = "viewMSAQC.asp?docID=" + docID + "&rev=" + revision + "&approvalType=" + approvalType;
    if (refreshPage)
		name += "&refresh=1";
	else
		name += "&refresh=0";
    openCenterWin(name,"viewMSAQC",1,1,0,0,0,0);
}

function popupRCMSFileResults(docID, revision, fileID, msid, approvalType, isCompanionFile)
{
    var url = "viewRCResults.aspx?docID=" + docID + "&rev=" + revision + "&fileID=" + fileID + "&msid=" + msid;

    if (approvalType > 0)
		url += "&type=" + approvalType ;

	if (typeof(isCompanionFile) != "undefined")
		url += "&cf=" + isCompanionFile;
	else
		url += "&cf=0";

    openCenterWin(url,"viewMSRC",1,1,0,0,0,0);
}

function popupRCMSResults(docID, revision, msid, approvalType)
{
    var url = "viewRCResults.aspx?docID=" + docID + "&rev=" + revision + "&msid=" + msid;

    if (approvalType > 0)
		url += "&type=" + approvalType ;

    openCenterWin(url,"viewMSRC",1,1,0,0,0,0);
}

function popupRCAmendedFiles(fileID, guid)
{
	var name = "viewRCAmendedFiles.asp?fileID=" + fileID + "&guid=" + guid ;
    var aNewWindow = window.open
    (
        name, "view_amended_files",
        "resizable=yes,scrollbars=no,height=200,width=250"
    );
    aNewWindow.focus();
}

var vlid = 0;
var vlprod = false;

function popupLetter(letterID, letterSecurityID, prod, sgTask)
{
    var name = "viewLetter.aspx?id=" + letterID + "&lsid=" + letterSecurityID;
    vlid = letterID;

    if (typeof(prod) != "undefined")
    {
    	name += "&prod=" + prod;
    	vlprod = prod;
    }
    else
    {
		vlprod = false;
    }

    if(typeof(sgTask) != "undefined")
    {
        name += "&sgTask=" + sgTask;
        sgTask = sgTask;
    }
    else
    {
        sgTask = false;
    }

    openCenterWin(name,"ViewLetter",1,1,0,0,0,0);
}


function popup_reviewer_details(docID, revision, ms_num, peopleID)
{
    var name = "reviewer_details.asp?docid=" + docID + "&rev=" + revision + "&ms_num=" + ms_num +
               "&peopleid=" + peopleID;
    openCenterWin(name,"reviewer_details",1,1,0,0,0,0);
}

function viewPrevVerPDFs(docID, msNum, revision, msid)
{
    if (revision == 0)
    {
    	var filename = msNum + ".pdf" ;
        var url = "fileDownload.asp?msid=" + msid + "&docID=" + docID + "&file=" + filename + "&rev=" + revision ;
        top.location.href = url;
    }
    else
    {
        var name = "viewPrevVerPDFs.asp?sub_num=" + docID + "&ms_num=" + msNum + "&rev=" + revision ;
        openCenterWin(name,"view_prev_ver_pdfs",1,1,0,0,0,0);
    }
}

function popupPublishInfo(docID, ms_num)
{
		var name = "med_publish_info.asp?docID=" + docID + "&ms_num=" + ms_num;
		openCenterWin(name,"publish_information",1,1,0,0,0,0);
}


function redirectToSubmitRevision(queryString)
{
	document.location.href = "ChooseEditSubmissionMethod.aspx?" + queryString;
}

function undoPropEditorAssignment(docID, revision, nChain, nAssignedBy, nEditorUnassigned)
{
    var name = "";
    name += "undoEditorAssignment.asp?docid=" + docID + "&rev=" + revision + "&chain=" + nChain;
    name += "&assignedBy=" + nAssignedBy  + "&editorUnassigned=" + nEditorUnassigned;

	openCenterWin(name,"undoEditorAssignment",1,1,0,0,0,0);
}

function undoEditorAssignment(docID, sessionThreadId)
{
	openCenterWin3("unassignEditor.aspx?documentid=" + docID + "&SessionThreadIdField=" + sessionThreadId ,"undoEditorAssignment",1,1,0,0,0,0);
}

function viewHelp(helpPage)
{
    openCenterWin(helpPage,"help",1,1,0,0,0,0);
}

function viewTransmittalForm(docID, revision, ms_num)
{
		var name = "";
		name += "transmittalForm.asp?docid=" + docID + "&rev=" + revision + "&ms_num=" + ms_num;
		name += "&includeMetadata=1" + "&isPopUpTransmittalForm=True";
		openCenterWin(name,"transmittalForm",1,1,0,0,0,0);
}


    function ViewOtherAuthorStatusPopUp(url)
    {
		var aNewWindow = window.open(url, "ViewOtherAuthorStatus","resizable=yes,scrollbars=yes,height=500,width=1010");
		
		var width = 1010;
		var height = 500;

		var x = (screen.availWidth - width) / 2; 
		var y = (screen.availHeight - height) / 2 

		aNewWindow.moveTo(x,y);
		
		aNewWindow.focus();
    }
   
    

    function ViewAuthorCorrHistoryPopUp(url)
    {
		var aNewWindow = window.open(url, "ViewAuthorCorrHistory","resizable=yes,scrollbars=yes,height=500,width=725");
		
		var width = 725;
		var height = 500;

		var x = (screen.availWidth - width) / 2; 
		var y = (screen.avaiHeight - height) / 2 

		aNewWindow.moveTo(x,y);
		
		aNewWindow.focus();
    }
   
    


function viewFieldHelp(theLink, windowName)
{
	openCenterWin(theLink,windowName,1,1,0,0,0,0);
}

function viewTransmittalHistory(docID, ms_num)
{
    var name = "";
    name += "transmittalFormHistory.asp?docid=" + docID + "&ms_num=" + ms_num;

	openCenterWin(name,"transmittalFormHistory",1,1,0,0,0,0);
}

function viewSavedTransmittalForm(transmittalHistoryID, ms_num)
{
    var name = "";
    name += "transmittalFormHistoryDetail.asp?";
    name += "transmittalHistoryID=" + transmittalHistoryID + "&ms_num=" + ms_num;

	openCenterWin(name,"transmittalFormHistoryDetail",1,1,0,0,0,0);

}

function viewTransmittalFormMetadata(docID, revision, manuscriptNumber)
{
    var name = "";
    name += "tfDownloadMetadata.asp?";
    name += "docid=" + docID + "&rev=" + revision + "&ms_num=" + manuscriptNumber;

	openCenterWin(name,"tfDownloadMetadata",1,1,0,0,0,0);
}

function checkAll(theForm, checkName, state)
{
    var checkCollection = eval("document.forms." + theForm + "." + checkName);

    if (typeof(checkCollection.length) == "undefined")
        checkCollection.checked = state;
    else
    {
        for (var k = 0; k < checkCollection.length; k++)
            checkCollection[k].checked = state;
    }
}

function showInfo(peopleID, docID, path, jrnlID)
{
    var url = path + "reviewer_info.asp?id=" + peopleID + "&docid=" + docID;

	if (typeof(jrnlID) != "undefined")
		url += "&jrnlID=" + jrnlID;

	openCenterWin(url,"peopleInfo",1,1,0,0,0,0)
}



function selectAll(theForm, selectName, state)
{
    var selectObj = eval("document.forms." + theForm + "." + selectName + ".options");

    if (typeof(selectObj.length) == "undefined")
        selectObj.selected = state;
    else
    {
        for (var k=0; k < selectObj.length; k++)
            selectObj[k].selected = state;
    }
}

function escapeDoubleQuotes(sUnescaped)
{
    var sEscaped = "" + sUnescaped.substring(0, sUnescaped.length);

    while (sEscaped.indexOf("\"") != -1)
        sEscaped = sEscaped.replace("\"", "&quot;");

    return sEscaped;
}

function unEscapeDoubleQuotes(sEscaped)
{
    var sUnescaped = "" + sEscaped.substring(0, sEscaped.length);

    while (temp.indexOf("&quot;") != -1)
        sUnescaped = sUnescaped.replace("&quot;", "\"");

    return sUnescaped;
}


function isNum(str)
{
	
	var strIdx = 0;

	if (str.charAt(0) == "-")
		strIdx++;

	while (strIdx < str.length)
	{
		if (!isDigit(str.charAt(strIdx)))
			return false;

		strIdx++;
	}

	return true;
}

function viewReviewerDisplayPreferences()
{
    var url;
    url = "reviewerDisplayPrefs.asp";

	openCenterWin(url,"viewReviewerDisplayPreferences",1,1,0,0,0,0);
}



function peopleURLpopup(theURL)
{
	openCenterWin(theURL,"peopleURLpopup",1,1,0,0,0,0);
}


function viewPersonalClassHelp()
{
    var url;
    url = "reviewerPersonalClassHelp.asp";

	openCenterWin(url,"viewPersonalClassHelp",1,1,0,0,0,0);
}

function viewSubmissionFiles(docID, revision, ms_num, approvalType)
{
    var url = "";
    url += "downloadSubFilesList.asp";
    url += "?docID=" + docID;
    url += "&docRevision=" + revision;
    url += "&msNum=" + ms_num;
    url += "&approvalType=" + approvalType;

    
	clientOpenCenterWinPercentage(url,"viewSubmissionFiles",1,1,0,0,0,0, 750, 500, 0.85);   
}



function viewAuthorInvLetter(authorInvLetterCorrHistID, letterSecurityID)
{
	vlid = authorInvLetterCorrHistID;
	var name = "viewLetter.asp?id=" + authorInvLetterCorrHistID + "&lsid=" + letterSecurityID;
	openCenterWin(name,"viewLetter",1,1,0,0,0,0);
}


function viewInvitationInfo(documentID, invitedId, idType)
{
    
	var name = "dotnetpopups/addtnlInvitationDetails.aspx" ;
	name += "?docID=" ;
	name += documentID ;
	name += "&invitedId=" ;
	name += invitedId ;
	name += "&idType=" ;
	name += idType ;
	openCenterWin(name,"ViewInvitationDetails",1,1,0,0,0,0);
}


function viewProposalDetails(docID, ms_num, sectionID)
{
    var intPos;

    var pageToLookFor = "editorCommentariesPending";

    var referringPage = document.location.href;

    intPos = referringPage.indexOf(pageToLookFor);

    if (intPos != -1)
    {
       sectionID = 1
    }

    var name = "";

    name += "EMDetails.aspx?";
    name += "docid=" + docID + "&ms_num=" + ms_num + "&sectionID=" + sectionID;

	openCenterWin(name,"proposal_details",1,1,0,0,0,0);
}


function doPreviewLetter(theURL)
{
	openCenterWin(theURL,"preview_letter",1,1,0,0,0,0);
}

function openAdditionalManuscriptDetails(sUrl, sName, sFeatures, bReplace)
{
	    var Win = window.open(sUrl,sName, sFeatures, bReplace);
	    Win.focus();
}



function openWinGeneric(sUrl, sName, sFeatures, bReplace)
{
	    var Win = window.open(sUrl,sName, sFeatures, bReplace);
	    Win.focus();
	  
} 


function OpenWindowFixedSize(url,name,width,height, scrollbar)     
{   
     var options = "width=" + width;   
     options += ",height=" + height;  
     options += ",scrollbars=" + scrollbar;        
     options += ",location=no,toolbar=no, menubar=no, resizable=no";   
  
     window.open(url,name,options);   
}   

function openCenterWin(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	openCenterWinPct(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.75);
}


function openCenterWinRet(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	return openCenterWinPctRet(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.75);
}


function openCenterWinPct(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct)
{
    openCenterWinPctRet(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct);
}


function openCenterWinPctRet(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct, addWidthHeightToUrl)
{

	var scrollArg = (scroll == false) ? '' : ',scrollbars=1';
	var resizeArg = (resize == false) ? '' : ',resizable=1';
	var menubarArg = (menubar == false) ? '' : ',menubar=1';
	var locationArg = (location == false) ? '' : ',location=1';
	var statusArg = (status == false) ? '' : ',status=1';
	var toolbarArg = (toolbar == false) ? '' : ',toolbar=1';
	var minW = minWidth; 
	var minH = minHeight; 

	var minW_hit = false; 
	var minH_hit = false; 
	var winW = 0; 
	var winH = 0; 
	var winL = 0; 
	var winT = 0; 
	var availW = 0; 
	var availH = 0; 
	if (navigator.userAgent.indexOf('Opera') != -1)
	{
		var minW = minWidth * 0.9; 
		var minH = minHeight * 0.8; 
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

		if (minW > availW - 70) {minW = availW - 70}
		if (minH > availH - 70) {minH = availH - 70}
	}

	else
	{
		availW = screen.availWidth; 
		availH = screen.availHeight; 
	}

	winW = availW * pct; 
	winH = availH * pct; 

	if (winW < minW) {winW = minW; minW_hit = true;}
	if (winH < minH) {winH = minH; minH_hit = true;}

	if (minW_hit) 
		winL = (availW - minW)/2;
	else
		winL = (availW - winW)/2;

	if (minH_hit) 
		winT = (availH - minH)/2;
	else
		winT = (availH - winH)/2;

    winW = Math.round(winW);
    winH = Math.round(winH);
    if (addWidthHeightToUrl && url != "")
    {
        url = url + "&width=" + winW + "&height=" + winH;
    }
	var Win = window.open(url,winName,"width="+winW+",height="+winH+scrollArg+resizeArg+menubarArg+locationArg+statusArg+toolbarArg);
	Win.focus();
	return Win; 
}


function openCenterWinStandardRet(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar)
{
	openCenterWin2(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar)
}


function openCenterWin2(url, winName, winW, winH, scroll, resize, menubar, location, status, toolbar)
{

	var scrollArg = (scroll == false) ? '' : ',scrollbars=1';
	var resizeArg = (resize == false) ? '' : ',resizable=1';
	var menubarArg = (menubar == false) ? '' : ',menubar=1';
	var locationArg = (location == false) ? '' : ',location=1';
	var statusArg = (status == false) ? '' : ',status=1';
	var toolbarArg = (toolbar == false) ? '' : ',toolbar=1';

	var winL = 0; 
	var winT = 0; 
	var availW = 0; 
	var availH = 0; 
	if (navigator.userAgent.indexOf('Opera') != -1)
	{

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

	else
	{
		availW = screen.availWidth; 
		availH = screen.availHeight; 
	}

	winL = (availW - winW)/2;

	winT = (availH - winH)/2;

	var Win = window.open(url,winName,"width="+winW+",height="+winH+",top="+winT+",left="+winL+scrollArg+resizeArg+menubarArg+locationArg+statusArg+toolbarArg);
	Win.focus();

	return Win;
}


function openCenterWin3(url, winName, scroll, resize, menubar, location, status, toolbar)
{
	openCenterWinPct(url, winName, scroll, resize, menubar, location, status, toolbar,
					 750, 500, 0.90);
}

function unauthorizedHTMLCheck(str)
{
	var scriptRegExp = new RegExp("<(/)*script.*>", "ig");
	if (str.match(scriptRegExp) != null)
		return("SCRIPT");

	var objectRegExp = new RegExp("<(/)*object.*>", "ig");
	if (str.match(objectRegExp) != null)
		return("OBJECT");

	var embedRegExp = new RegExp("<(/)*embed.*>", "ig");
	if (str.match(embedRegExp) != null)
		return("EMBED");

	var appletRegExp = new RegExp("<(/)*applet.*>", "ig");
	if (str.match(appletRegExp) != null)
		return("APPLET");

	return("");
}


function trimCR(str)
{
	str = escape(str)

	var i = -1;

	if (str.indexOf("%0D%0A") != -1)
	{
		i = str.indexOf("%0D%0A");
		while(i != -1)
		{
			str = str.substring(0, i) + " " + str.substring(i+6, str.length);
			i = str.indexOf("%0D%0A");
		}
	}

	else if (str.indexOf("%0A") != -1)
	{
		i = str.indexOf("%0A");
		while(i != -1)
		{
			str = str.substring(0, i) + " " + str.substring(i+3, str.length);
			i = str.indexOf("%0A");
		}
	}

	else if (str.indexOf("%0D") != -1)
	{
		i = str.indexOf("%0D");
		while(i != -1)
		{
			str = str.substring(0, i) + " " + str.substring(i+3, str.length);
			i = str.indexOf("%0A");
		}
	}
	return unescape(str);
}



function disableAllButtons(theForm)
{
	for (var i = 0; i < theForm.elements.length; i++)
	{
		if (theForm.elements[i].type == 'button' || theForm.elements[i].type == 'submit')
			theForm.elements[i].disabled = true;
	}
}


function setLinkState(state, menu, content)
{
	var a = 0
	if (menu)
	{
		if(window.top.menuPage)
		{
			while (window.top.menuPage.document.links[a])
			{
				window.top.menuPage.document.links[a].disabled = state;
				a++;
			}
		}
	}
	if (content)
	{
		a = 0;
		while (window.document.links[a])
		{
			window.document.links[a].disabled = state;
			a++;
		}
	}
}


function checkLinkStateWithNoAlerts(link)
{
	if (link.disabled == true)
		return false;
	else
		return true;
}


function trimAllTextFields(theForm)
{
	for (var i = 0; i < theForm.elements.length; i++)
	{
		if (theForm.elements[i].type == 'text' || theForm.elements[i].type == 'textarea')
			theForm.elements[i].value = trim(theForm.elements[i].value);
	}
}


function areEqual(variableOne, variableTwo)
{
	return (variableOne == variableTwo ? true : false );
}


function isEmptyString(theString)
{
	return  ( (trim(theString)) == "" ? true : false );
}


function doesObjectExist(theObject)
{
	return ( theObject != null );
}


function popupScheduleGroupDetails(scheduleGroupID)
{
	var name = "";
    name += "ScheduleGroupDetails.aspx?";
    name += "scheduleGroupID=" + scheduleGroupID;
    
    openCenterWin(name,"popupScheduleGroupDetailsWindow",1,1,0,0,0,0);
}


function parseBoolean(theValue)
{
	var theString = ("" + theValue).toLowerCase();

	if(theString == "true" || theString == "1")
		return true;

	return false;
}


function isSomethingChecked(theForm)
{
	for (var i = 0; i < theForm.elements.length; i++)
	{
		if (theForm.elements[i].type == 'checkbox')
		{
			if( parseBoolean(theForm.elements[i].checked) )
				return true;
		}
	}

    return false;
}


function checkOffAll(state)
{
	$('input').prop('checked', state);
}


function popupTechCheckWindow(docID)
{
    var name = "technicalCheck.aspx?docid=" + docID;
    openCenterWinPct(name,"technicalCheck",1,1,0,0,0,0, 600, 600, 0.75);
}


function popupTechComments(techHistID)
{
    var name = "viewTechComments.asp?id=" + techHistID;
    openCenterWinPct(name,"technicalComments",1,1,0,0,0,0, 600, 600, 0.75);
}


function popupAuthorsReviewerPreferences(docID, revision)
{
    var name = "viewAuthorsReviewerPreferences.aspx?docid=" + docID + "&revision=" + revision;
    openCenterWin(name,"popupAuthorsReviewerPreferencesWindow",1,1,0,0,0,0)
}

function popupAuthorResponseToReviewers(docID, revision)
{
    var name = "viewAuthorsResponseToReviewers.aspx?docid=" + docID + "&revision=" + revision;
    openCenterWin(name,"popupAuthorResponseToReviewersWindow",1,1,0,0,0,0)
}

function isAlpha(ch)
{
	return ((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z'));
}

function isDigit(ch)
{
	return(!isNaN(parseInt(ch)));
}


function verifyEmailAddresses(addresses)
{
	var addrs = addresses.split(";");
	
	for(var x = 0; x < addrs.length; x++)
	{
		var trimmedStr = trim(addrs[x]);
		if(!verifyEmailAddr(trimmedStr))
		{
			return false;
		}
	}
	
	return true;
}

function popupEditItemMetadata(id, working, filetypeid, description)
{
    var name = "editSubmissionItemMetadata.aspx?tform=0&id=" + id + (working == true ? "&working=1" : "") ;
    if (filetypeid !== undefined)
    {
        name += "&filetypeid=" + filetypeid;
    }
    if (description !== undefined)
    {
        name += "&description=" + encodeURIComponent(description);
    }
	openCenterWin2(name, "name", 800, 400, 1, 1, 0, 0, 0, 0);
}


function popupChangeColorForFile(companionFileID, fileGroupName, fileName, parentHTMLColumnID)
{
	var changeColorURL = "changeColorCode.asp?" +
						 "companionFileID=" + companionFileID +
						 "&fileTypeName=" + fileGroupName +     
						 "&fileName=" +fileName +
						 "&parentHTMLColumnID=" + parentHTMLColumnID;
	openCenterWin(changeColorURL, "changeColorWindow", 1, 1, 0, 0, 0, 0);
    return false;
}


function isValidColorCodeValue(colorCodeValue)
{
	
	return ( colorCodeValue.match(/^\#(\d|[A-F]|[a-f]){6}$/) );
}


function validateThisDate(inputHoldingDate, inputHoldingDateCanBeEmpty, alertMessageIfEmpty, alertMessageIfDateInPast)
{
	var isValidDate = true;

	inputHoldingDate.value = trim(inputHoldingDate.value);

	if( (inputHoldingDate.value == "") && !inputHoldingDateCanBeEmpty)
	{
		alert(alertMessageIfEmpty);
		return false;
	}
	else if(!convert_date(inputHoldingDate) )
	{
		return false;
	}
	else if(
				(inputHoldingDate.value != "") &&
				( new Date(inputHoldingDate.value) < (new Date()).setHours(0,0,0,0) )
			)
	{
		alert(alertMessageIfDateInPast);
		return false;
	}

	return true;
}

	function addWorkingDays(startDate, addDays)
	{
		startDate = new Date(startDate);
		var n = 0;
		var offset = addDays;
		if (offset > 0)
		{
			for (i=1; n<offset; i++)
			{
				var tmp = new Date(startDate.getFullYear(),startDate.getMonth(),startDate.getDate()+i);
				if (tmp.getDay() != 0 && tmp.getDay() != 6){n++}
			}
			i--;
		}
		else
		{
			for (i=1; n>offset; i++)
			{
				var tmp = new Date(startDate.getFullYear(),startDate.getMonth(),startDate.getDate()-i);
				if (tmp.getDay() != 0 && tmp.getDay() != 6){n--}
			}
			i = (i*-1)+1;
		}
		var resultDate = new Date(startDate.getFullYear(),startDate.getMonth(),startDate.getDate()+i);
		var slashDate = resultDate.getMonth()+1+"/"+resultDate.getDate()+"/"+resultDate.getFullYear();
		return(slashDate);
	}


	function fixIESelectionBug()
	{
		if (window.createPopup && document.compatMode && document.compatMode=="CSS1Compat")
		{
			document.onreadystatechange = onresize = function fixIE6AbsPos()
			{
				if (!document.body) return;
				if (document.body.style.margin != "0px") document.body.style.margin = 0;
				onresize = null;
				document.body.style.height = 0;
				setTimeout(function(){ document.body.style.height = document.documentElement.scrollHeight+"px"; }, 1);
				setTimeout(function(){ onresize = fixIE6AbsPos; }, 100);
			}
		}
	}


function popupSetFinalDisposition(docID, revision, ms_num)
{

    var name = "setFinalDisposition.aspx?docid=" + docID + "&rev=" + revision + "&ms_num=" + ms_num;
    openCenterWin(name,"setFinalDispotion",1,1,0,0,0,0);
}

function popupClassificationDetail(classID, classJournalCode)
{
    var name = "DotNetPopUps/ViewClassificationDetails.aspx?classcid=" + classID;
    
    if (typeof(classJournalCode) != "undefined")
    {
    	name += "&classJournalCode=" + classJournalCode;
    }
    
    openCenterWin2(name, "ViewClassificationDetails" + classID, 850, 200, 1, 1, 0, 0, 0, 0);
}

function popupEditPersonalKeyword(keywordPeopleDataID, keywordID, keywordIDList)
{
    var name = "DotNetPopUps/editPersonalKeyword.aspx?keywordPeopleDataID=" + keywordPeopleDataID + "&keywordID=" + keywordID;
    if (typeof(keywordIDList) != "undefined")
    {
    	name += "&keywordIDList=" + keywordIDList;
    }
    openCenterWin2(name, "EditPersonalKeyword", 850, 200, 1, 1, 0, 0, 0, 0);
}


function OpenPopUpGeneral(sUrl, sWinName, sFeatures, bReplace)
{
	var winOpen = window.open(sUrl, sWinName, sFeatures, bReplace);
	winOpen.focus();
}

function cleanEmailString(emailString)
{

	if (typeof emailString != 'string')
	{
		return "";
	}

	var re = /^[\s;]*/;
	var cleanedEmailString = emailString.replace(re, "");

	re = /[;\s]*$/;
	cleanedEmailString = cleanedEmailString.replace(re, "");


	re = /[;\s]+/g;
	cleanedEmailString = cleanedEmailString.replace(re, ";");
	
	return cleanedEmailString;
}


function ToggleCollapsedText(divToExpCollapse, spanClicked, expandedText, collapsedText)
{
	var spanElement = document.getElementById(spanClicked);
	var elementToExpCollapse = document.getElementById(divToExpCollapse).style; //get the elementToExpCollapse
	
	if(elementToExpCollapse.display == "none")
	{
		elementToExpCollapse.display = "";
		spanElement.innerHTML = expandedText;
	} 
	else 
	{
		elementToExpCollapse.display = "none";
		spanElement.innerHTML = collapsedText;
	}
}

function WarnRedirect(message, link)
{
	if (window.confirm(message))
	document.location.href = link;
}
function redirectWarning(message, link)
{
    if (typeof timerID != 'undefined') clearTimeout(timerID);
    var okAction = function()
    {
        document.location.href = link;
    };
    var cancelAction = function()
    {
        resetTimeout();
    };
    $("#warningDialog")[0].showDialog(message, okAction);
}
	
function SetFormElementValue(formElementID, sValue)
{
	var elementToSet = document.getElementById(formElementID);
	
	if(elementToSet != null)
	{
		elementToSet.value = sValue;
	}			
}


function getScheduleGroupInventoryPopUp(scheduleGroupID, fromPage)
{
    var url = "";
    url += "ScheduleGroupFileInventory.aspx";
    url += "?scheduleGroupID=" + scheduleGroupID;
    url += "&fromPage=" + fromPage;
    
	openCenterWin(url,"viewSGFileInventory",1,1,0,0,0,0);
}

function BrowserIsIE()
{
	var isIE;		// indicates whether Internet Explorer is being used
	if (navigator.appName == "Microsoft Internet Explorer" || !!navigator.userAgent.match(/Trident\/7\./))  
	{
		isIE = true;
	}
	else
	{
		isIE = false;
	}
	
	return isIE;
}



Date.prototype.toStandardString = function()
{
	
	if (!this.valueOf())
	{
		return this;
	}
	
	var dt = this;
	
	
	var mm = ((dt.getMonth() < 9) ? "0" : "") + (dt.getMonth() + 1);
	
	
	var dd = ((dt.getDate() < 10) ? "0" : "") + dt.getDate();
	
	
	var yyyy = (dt.getFullYear() < 2000) ? (dt.getFullYear() + 100) : dt.getFullYear();
	
	
	var dateString = mm + "/" + dd + "/" + yyyy;
	
	return dateString;
}


Date.prototype.toSortableString = function()
{
	
	if (!this.valueOf())
	{
		return this;
	}
	
	var dt = this;
	
	var paddedYear = '' + (10000 + dt.getFullYear());
	var paddedMonth = '' + (100 + dt.getMonth() + 1);
	var paddedDate = '' + (100 + dt.getDate());
	
	return paddedYear.substring(1) + '-' + paddedMonth.substring(1) + '-' + paddedDate.substring(1);
}


function cleanUpDate(fieldID)
{
	var txtDate = document.getElementById(fieldID);
	var dateStr = txtDate.value;
	
	dateStr = trim(dateStr);
	dateStr = dateStr.replace(/^(\d)\//, '0$1/')		// Pad month
	dateStr = dateStr.replace(/\/(\d)\//, '/0$1/')		// Pad day
	dateStr = dateStr.replace(/\/(\d\d)$/, '/20$1')		// Pad year
	
	txtDate.value = dateStr;
}



function disableEnterKeySubmit(event)
{
   var e = (event ? event : window.event); 
   return (e.keyCode ? e.keyCode : e.charCode ? e.charCode : e.which) != 13;
}

function popupFeeList(docID)
{
	openCenterWin3("viewFeeList.aspx?docID=" + docID,"viewFeeList",1,1,0,0,0,0);	
}

function emailAddrListObj(emailAddrList)
{
    this.emailAddrList = emailAddrList;
    this.invalidList = "";
    this.numInvalid = 0;
}


function IsValidEmailAddressList(emailAddrObj)
{
    var invalidEmailAddr = new Array();
	var emailAddrArray = emailAddrObj.emailAddrList.split(";");
    

    for (var i = 0; i < emailAddrArray.length; i++)
    {
		
		var sEmailAddressAll = trim(emailAddrArray[i]);
		var emailAddrLessName = "";
		var emailAddrToValidate = "";
		
		if(sEmailAddressAll.length == 0)
		{
			continue;
		}

		var emailAddrRegExp = new RegExp(/^"[^"\r\n]*"/);

		if(emailAddrRegExp.test(sEmailAddressAll))
		{
			emailAddrLessName = sEmailAddressAll.replace(emailAddrRegExp, "");
		}
		

		emailAddrToValidate = ((emailAddrLessName.length != 0) ? emailAddrLessName : sEmailAddressAll);
		
		if(! IsValidEmailAddress(emailAddrToValidate))
		{
			invalidEmailAddr.push(emailAddrToValidate);
            emailAddrObj.numInvalid++;
		}
	}
	
    if (emailAddrObj.numInvalid > 0)
    {
        emailAddrObj.invalidList = invalidEmailAddr.join("; ");
        return false;
    }
	return true;	
}



function IsValidEmailAddress(emailToValidate)
{
	emailToValidate = trim(emailToValidate);
	var emailAddressRegexString = "[-a-zA-Z0-9!#$%&'*+/=?^_'{|}~]+(?:\\.[-a-zA-Z0-9!#$%&'*+/=?^_'{|}~]+)*" +		// Local part
				                   "@([a-zA-Z0-9](?:[-a-zA-Z0-9]*[a-zA-Z0-9])?(?:\\.[a-zA-Z0-9](?:[-a-zA-Z0-9]*[a-zA-Z0-9])?)+|" +	// Domain name or
				                   "\\[(?:\\d{1,3}(?:\\.\\d{1,3}){3}|IPv6:[0-9A-Fa-f:]{4,39})\\])";                 // IP Address
					                   
	var emailAddrRegExp = new RegExp(emailAddressRegexString);
	
	return (emailAddrRegExp.test(emailToValidate));
}


function isValidClassNumber(inputString) {

	var regExp = new RegExp("^(?:\.?[0-9]+)+$");

	return (inputString.match(regExp) == null);
}


function showToolTip(ctrl,text) 
{
	ctrl.title=text;
}


function clientOpenCenterWinPercentage(url, winName, scroll, resize, menubar, location, status, toolbar,
						  minWidth, minHeight, pct)
{
	// Previously "openCenterWinPct" in clientutil.js
	// copied from OpenPopWindowMain.js and renamed to not create conflict
	
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
	if (minW_hit) //If popup height will be less than minimum then use the minimum to center
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
	var Win = window.open(url,winName,"width="+winW+",height="+winH+",top="+winT+",left="+winL+scrollArg+resizeArg+menubarArg+locationArg+statusArg+toolbarArg);
	Win.focus();
	return Win;
}


function hasCssClass(element, className) {
    var rxp = new RegExp('\\b' + className + '\\b');
    return rxp.test(element.className);
}

function addCssClass(element, className) {
    if (!hasCssClass(element, className)) {
        if (element.className == '') {
            element.className = className;
        } else {
            element.className += ' ' + className;
        }
    }
}

function removeCssClass(element, className) {
    if (hasCssClass(element, className)) {
        var tmp = element.className;
        tmp = tmp.replace(className, '');
        tmp = tmp.replace(/\s{2,}/, ' ');
        tmp = tmp.replace(/^\s/, '');
        tmp = tmp.replace(/\s$/, '');
        element.className = tmp;
    }
}


function AddPopUpToCloseList(obj)
{
	var identWindow = null;
	var root = null;	        
	GetOpener(obj);

	function GetOpener(obj)
	{
	    if (obj.opener)
	    {
	        root = obj.opener;
	        GetOpener(root);
	    }
	}

	// Now we need to get the ident frame that holds the openWindows array
	if (root && root.top && root.top.ident)
	{
	    identWindow = root.top.ident;

	    if (identWindow.openWindows)
	    {
	        identWindow.openWindows[identWindow.openWindows.length] = window;
	    }

	}
}

function popupMedlineSearch(searchString) {
    var url = "https://www.kfinder.com/member-search/login.cgi?medweb=1&data=qhjnK2a9jJgT28s2GQY8YGwvX8XUOvW8W6pvj85npuq8hq&searchstring=";
    url += searchString;
    url += "&dbproduct=MEDLINE&searchlogic=fuzzy&getcount=200&relevance=50&segments=4&getchunk=20&concept_mapping=on&wordvars=on&relevance_sort=on";
    window.open(url);
}

function showPeopleInfo(id, docID)
{	
    var url = "reviewer_info.asp?id=" + id + "&docid=" + docID;
    clientOpenCenterWinPercentage(url, "showPeopleInfo", 1, 1, 0, 0, 0, 0,750, 500, 0.75);
}

function ShowElement(elmName, show)
{
    document.getElementById(elmName).style.display = show ? 'block' : 'none';
    return false;
}


var numPosts = 0;
function PreventDoublePost() 
{
    if (document.activeElement.tagName.toLowerCase() !== 'button')
    {
        return true;           
    }
    else if(numPosts === 0)
    {
        numPosts++;
        return true;
    } 
    else 
    {
        return false;
    }
}
