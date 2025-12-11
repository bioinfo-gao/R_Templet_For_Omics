
function popupViewRevAttach(docID,reviewerAttachmentViewContext)
{
    // 11.0-28 20131108 JonC    Changed URL to point to newly ported dot net page Editor_ViewRevAttach.aspx. 
	var url = "Editor_ViewRevAttach.aspx?docid=" + docID + "&reviewerAttachmentViewContext=" + reviewerAttachmentViewContext;

    var aNewWindow = window.open(url,'Editor_ViewRevAttach','resizable=1,scrollbars=1,height=550,width=750');
  
    aNewWindow.focus();
}

function SetAllCheckboxes(theForm, value, szCheckboxFieldName)
{
	var iCheckbox = 0;
	var i = 0;

	for(i = 0; i < theForm.elements.length; i++)
	{
		if (theForm.elements[i].type == "checkbox" && theForm.elements[i].name == (szCheckboxFieldName + "_" + iCheckbox))
		{
			theForm.elements[i].checked = value;
			iCheckbox++;
		}
	}

}