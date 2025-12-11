

// Prompt for Remove the Reviewer Attachments
//	This version uses the "retForm" variable 
// #26835   20140701    JonC. Added the szparams variable to the parameter list for the Remove URL. 
function RemoveReviewerAttachment(attachmentid, szparams)
{
    if (window.confirm("Are you sure you want to permanently remove this item\x3f"))
    {
    	document.retForm.action='Remove_RevAttach.ASP?attachmentid=' + attachmentid + szparams ;
		document.retForm.submit();
    }
}

// Prompt for Remove the Reviewer Attachments
//	This version uses the "revForm" variable
function RemoveReviewerAttachment2(attachmentid, szparams)
{
    if (window.confirm("Are you sure you want to permanently remove this item\x3f"))
    {
    	document.revForm.action='Remove_RevAttach.ASP?attachmentid=' + attachmentid;
		document.revForm.submit();
    }
}



