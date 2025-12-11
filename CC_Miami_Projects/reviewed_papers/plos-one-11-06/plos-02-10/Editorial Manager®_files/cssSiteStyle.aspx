

/**********************************************************************************************************
********
********  GENERAL STYLES
********
***********************************************************************************************************/
/*****************************************************************
	Styles for main document and 2-column containers 
******************************************************************/
body
{
	width: 97%; /*Set to 97% to fix horizontal scroll bug in 800x600*/
	padding-right: 0px;
	padding-left: 7px;
	padding-top: 20px;
	padding-bottom:20px;
	margin: 0px;
	background-color: #ffffff;
	color: #000033;
	font-size: 70%;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	clear: both;
}
body.admin
{
    
}
body.report
{
    
}
body.implMessageBanner
{
	padding: 0px;
	margin: 0px;
	background-color: #ff0000;
	color: #ffffff;
	font-weight: bold;
	font-size: 100%;
	font-family:arial,helvetica,sans-serif;
}
body.proxyBanner
{
	padding: 0px;
	margin: 0px;
	background-color: #fcee00;
	color: #000000;
	font-size: 100%;
	font-family:arial,helvetica,sans-serif;
}
body.emLogo
{
	padding: 0px;
	margin: 0px;
}
.emLogoBackground
{
	width: 100%;
	background-color: #CCCCCC;
	color: #000000;	
	font-size: 8pt;
	font-family:arial,helvetica,sans-serif;
}
div.left_side
{
    padding-right: 10px;
	padding-left: 15px;
	padding-bottom: 10px;
	padding-top: 10px;
	width: 220px;
	
    float: left;
}
div.left_side_absolute
{
	padding-right: 10px;
	padding-left: 15px;
	padding-bottom: 10px;
	padding-top: 10px;
	width: 180px;
    position: absolute; 
}
div.right_side
{
	padding-top: 0px;
	padding-right: 10px;
	padding-left: 20px;
	padding-bottom: 10px;
	border-left: #cccccc 1px dotted;
	width: 550px;
	float: left;

}
div.right_side_with_no_left_border
{
	padding-right: 10px;
	padding-left: 20px;
	left: 240px;
	padding-bottom: 10px;
	width: 550px;
	padding-top: 0px;
	position: absolute;
}
div.left_side_with_results_table
{
	padding-right: 10px;
	padding-left: 15px;
	padding-top: 10px;
	width: 220px;
	float: left;
}
div.right_side_with_results_table
{
	padding-top: 0px;
	padding-right: 10px;
	padding-left: 20px;
	padding-bottom: 10px;
	border-left: #cccccc 1px dotted;
	width: 550px;
	float: left;
}
div.left_side_nowrap 
{
	padding-right: 10px;
	padding-left: 15px;
	padding-bottom: 10px;
	padding-top: 10px;
	width: 180px;
	float:left;
}
div.right_side_nowrap
{
	padding-top: 0px;
	padding-right: 10px;
	padding-left: 20px;
	padding-bottom: 10px;
	border-left: #cccccc 1px dotted;
	width: 550px;
	margin:-3px 0px 0px -2px;
	position:relative;
	float:left;
}
div.background_div
{
	border: #cccccc 1px solid;
	padding-right: 10px;
	padding-left: 20px;
	padding-top: 10px;
	padding-bottom: 10px;
	background-color: #f4f4f4;
}
div.popupFormContainer
{
	padding-left: 20px;
	padding-right: 20px;
	width: 95%;
}
div.pageTitle, span.pageTitle, p.pageTitle
{
	font-size: 150%;
	font-weight: bold;
    text-align: center;
}
p.instructions, div.instructions
{
	text-align:left;
	margin-left:1.5cm;
	margin-right:1.5cm
}
div.topheader
{
	text-align: center;
	margin-top:10px;
	margin-bottom:-10px;
	padding-right: 10px;
	padding-left: 20px;
	width: 97.7%;
}
div.topheaderNoBottomMargin
{
	text-align: center;
	margin-top:10px;
	padding-right: 10px;
	padding-left: 20px;
	width: 97.7%;
}

/*****************************************************************
	Styles for fieldsets and legends 
*****************************************************************/
.fieldset_center
{
	width: 85%; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
	display: block;
}
.fieldset_invisible
{
	position: relative;
	margin-top: 0px;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 20px;
}
.fieldset_not_visible
{
	position: relative;
	margin-top: 20px;
	border: #f4f4f4 1px solid;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 10px;
}
.fieldset_nolegend
{
	position: relative;
	margin-top: 0px;
	border: #cccccc 1px solid;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 20px;
	background-color: #f4f4f4;
}
.fieldset_nolegend_lesspadding
{
	position: relative;
	margin-top: 0px;
	border: #cccccc 1px solid;
	padding: 2px;
	background-color: #f4f4f4;
}
fieldset, .divFieldset
{
	position: relative;
	margin-top: 0px;
	border: #cccccc 1px solid;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	background-color: #f4f4f4;
}
.fieldset2
{
	position: relative;
	margin-top: 0px;
	border: #cccccc 1px solid;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	background-color: #FAFAFA;
}
legend, .divFieldset > #legend
{
	position: relative;
	margin-top:-10px;
	margin-bottom:7px;
	background:#ffffff;
	border:1px solid #cccccc;
	padding: 2px 5px;
	font-weight: bold;
	text-align:left;
}
legend a
{
    color:#0000FF;
    text-decoration:underline;
    cursor:pointer;
}

/*****************************************************************
	Styles for page "containers"
*****************************************************************/
div.centered, table.centered
{
    margin-left: auto;
    margin-right: auto;
}
div.faqcontainer
{
	border: 1px solid #003365; 
	padding: 10px; 
	width: 550px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #ffffff;
}
div.faqbox
{
	padding: 10px;
	border-top: #002233 1px solid;
	border-bottom: #002233 1px solid;
	left: 10px;
	width: 490px;
	background-color: #ffffff;
}
div.boundedInstructions
{
	padding: 5px; 
	border-top: #002233 1px solid;
	border-bottom: #002233 1px solid;
	background-color: #ffffff;
}
div.refineReportInstructions
{
	margin-top:5px;
	margin-bottom:5px;
	margin-right:10px;
}
div.reviewerCommentsBox
{
	padding: 10px;
	border: #002233 1px solid;
	background-color: #ffffff;
}
div.centeredboxlesspadding
{ 
	padding: 2px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
}
div.centeredboxlesspaddingleft
{ 
	padding: 2px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
    background-color: #f4f4f4;
	text-align:left;
	line-height:14px;
}
div.centeredboxlesspaddingleft_nobackground
{ 
	padding: 2px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	text-align:left;
	line-height:14px;
}
div.centeredbox
{ 
	padding: 10px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
}
div.centeredboxalt
{ 
	padding: 10px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	line-height:14px;
}
div.centeredboxlarge
{ 
	padding: 10px; 
	width: 90%; 
	max-width: 780px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:center;
}
div.searchResultlarge
{ 
	padding: 10px; 
	width: 90%; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:center;
}
div.centeredboxlargenoborder
{ 
	padding: 10px; 
	width: 90%; 
	max-width: 780px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:center;
}
div.centeredboxlargewithborder
{ 
	padding: 0px; 
	width: 100%; 
	max-width: 1000px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:left;
    border: 1px solid #cccccc; 
}
div.centeredboxlargealt
{ 
	padding: 10px; 
	width: 90%; 
	max-width: 600px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:center;
}
div.scrollboxDiv
{ 
	padding: 2px; 
	width: 300px; 
	height: 400px;
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:left;
	overflow:auto;
}
div.ShadedLeftAlignDiv
{ 
	padding: 2px; 
	line-height:14px;
	text-align:left;
}
fieldset.centeredboxForControlFieldSet
{ 
	padding: 2px; 
	width: 650px; 
	margin-top: 10px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #CCCCCC;
	line-height:14px;
	text-align:center;
}
.panelShadedWide
{
	background-color: #CCCCCC;
	width: 780px;
	text-align:center;
}
div.centeredboxlargenowidth
{ 
	padding: 10px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:1.4em ;
	text-align:center;
}
div.centeredboxlargenowidthNoBorder
{ 
	padding: 10px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
    line-height:1.4em; 
	text-align:center;
}
div.centeredboxlargeLeft
{ 
	padding: 10px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:left;
}
div.centeredboxLeft
{ 
	padding: 10px; 
	width: 90%; 
	max-width:650px;
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	background-color: #f4f4f4;
	line-height:14px;
	text-align:left;
}
div.centeredboxlargenopadding
{ 
	width: 90%; 
	max-width: 750px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	line-height:14px;
}
div.centeredboxwidthofscreen
{
	margin-top: 0px;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 10px;
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	text-align: center;
}
div.enclosedWhiteBox
{
	padding: 10px; 
	margin-top: auto; 
	margin-left: auto; 
	margin-right: auto; 
	margin-bottom: auto; 
	border: 1px solid #cccccc; 
	background-color: white;
	line-height:14px;
}
div.messageBox
{ 
	padding: 2px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto;
    margin-bottom: 5px;
	border: 1px solid #cccccc; 
	background-color: #f4f4f4;
	line-height:14px;
}

.PageHeader
{ 
	padding: 0px; 
	width: 90%; 
	max-width: 650px; 
	margin-top: 5px; 
	margin-left: auto; 
	margin-right: auto; 
	text-align:center;
    display:block;
}

/*****************************************************************
Element Styles for browser consistency
******************************************************************/
a:link
{
	text-decoration: none;
    color: #0000ee;
}
a:visited
{
	text-decoration: none;
	color:#0000ee;
}
a:hover
{
	text-decoration: underline
}
input[type="text"], input[type="password"], input:not([type]), textarea, 
input[type="button"].txt_disabled, input[type="button"].req_disabled
{
    border-width: 1px;
    border-color: #f0f0f0;
    border-style: inset;
    border-radius: 0px;
    padding: 2px;
    margin: 0;
    font-size: 100%;
    font-family: verdana, geneva, arial, helvetica, sans-serif;
}
select
{
    border-width: 1px;
    border-color: #f0f0f0;
    border-style: inset;
    border-radius: 0;
    padding: 0px 0px 0px 1px;
    font-size: 95%;
    font-family: verdana, geneva, arial, helvetica, sans-serif;
	color: #000000;
	vertical-align: middle;
	margin-top: 2px;
	margin-bottom: 2px;
    cursor: pointer;
}
select:not([multiple]) 
{
    /* make arrow and background */
    background: linear-gradient(45deg, transparent 50%, #888 50%),
                linear-gradient(135deg, #888 50%, transparent 50%),
                linear-gradient(to right, #000, #e6e6e6 1px),
                linear-gradient(to left, #fff, #fff);
    background-position: calc(100% - 12px) calc(8px),
                         calc(100% - 8px) calc(8px),
                         100% 0,
                         100% 0;
    background-size: 4px 4px,
                     4px 4px,
                     24px 18px,
                     100% 18px;
    background-repeat: no-repeat;

    /* styling and reset */
    padding: 1px 25px 1px 1px;
    min-width: 36px;
    height: 18px;
    /* reset */
    -webkit-box-sizing: border-box;
    -moz-box-sizing: border-box;
    box-sizing: border-box;
    -webkit-appearance:none;
    -moz-appearance:none;
}
select::-ms-expand 
{
    display: none;
}
select:not([multiple]):hover 
{
    background: linear-gradient(45deg, transparent 50%, #454545 50%), 
                linear-gradient(135deg, #454545 50%, transparent 50%), 
                linear-gradient(to right, #000, #dadada 1px),
                linear-gradient(to left, #fff, #fff);
    background-position: calc(100% - 12px) calc(8px),
                         calc(100% - 8px) calc(8px),
                         100% 0,
                         100% 0;
    background-size: 4px 4px,
                     4px 4px,
                     24px 18px,
                     100% 18px;
    background-repeat: no-repeat;
}
input[type="button"], input[type="submit"], input[type="file"], button[type="button"]
{
    padding: 1px 5px;
    color: #003365;
	font-weight: normal;
    font-size: 100%;
    font-family: verdana, geneva, arial, helvetica, sans-serif;
    background-color: #cecece;
    background-image: linear-gradient(#e4e4e4 50%, #cecece 100%);
    border: 1px solid #cecece;
    border-radius: 0.35em;
    box-shadow: 1px 1px 2px #555555;
    cursor: pointer;
}
input[type="file"]
{
    border-top-left-radius: 0;
    border-bottom-left-radius: 0;
    outline: 0 !important;
}
input[type="text"].browse
{
    border-style: solid;
}
input[type="button"]:focus, input[type="file"]:focus, input[type="submit"]:focus, button[type="button"]:focus
{
    background-image: linear-gradient(#dedede 50%, #adadad 100%);
    border:1px solid #777777;
    outline: 0 !important;
}
input[type="button"]:hover, input[type="file"]:hover, input[type="submit"]:hover, button[type="button"]:hover
{
    background-image: linear-gradient(#dedede 50%, #adadad 100%);
    border:1px solid #777777;
    cursor: pointer;
}
input[type="button"]:active, input[type="file"]:active, input[type="submit"]:active, button[type="button"]:active
{
    background-color: #e4e4e4;
    background-image: none;
    box-shadow: inset 1px 1px 2px #555555;
    outline: 0 !important;
}
input[type="button"]:disabled, input[type="file"]:disabled, input[type="submit"]:disabled, button[type="button"]:disabled
{
    border:1px solid #cecece;
    box-shadow: 1px 1px 2px #777777;
    color:#999999;
}
input[type="button"]:hover:disabled, input[type="file"]:hover:disabled, input[type="submit"]:hover:disabled, button[type="button"]:hover:disabled
{
    border:1px solid #cecece;
    background-image: linear-gradient(#e4e4e4 50%, #cecece 100%);
    box-shadow: 1px 1px 2px #777777;
    color:#999999;
    outline: 0 !important;
}
/*****************************************************************
Styles and Style Classes for hyperlinks
*****************************************************************/
.linkInActiveXSmallBold
{
	font-size: x-small;
	color: #808080;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkInActiveXSmallBold:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkInActiveXSmallBold:active
{
	color: #9b0521;
}
.linkInActiveXSmall
{
	font-size: x-small;
	color: #808080;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkActive
{
	cursor: pointer;
	color:blue;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkInActiveXSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkInActiveXSmall:active
{
	color: #9b0521;
}
.linkBlackXSmall
{
	font-size: x-small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkBlackXSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkBlackXSmall:active
{
	color: #9b0521;
}
a.linkRedXSmall:link
{
	font-size: x-small;
	color: #ff0000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
a.linkRedXSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
a.linkRedXSmall:active
{
	color: #9b0521;
}
.linkBlackXXSmall
{
	font-size: xx-small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkBlackXXSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkBlackXXSmall:active
{
	color: #9b0521;
}
.linkBlackXSmallBold
{
	font-weight:bold;
	font-size: x-small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkBlackXSmallBold:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkBlackXSmallBold:active
{
	color: #9b0521;
}
a.linkRedXSmallBold:link
{
	font-weight:bold;
	font-size: x-small;
	color: #ff0000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
a.linkRedXSmallBold:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
a.linkRedXSmallBold:active
{
	color: #9b0521;
}
.linkBlackSmall
{
	font-size: inherit;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	vertical-align: middle;
}
.linkBlackSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
    text-decoration: none;
}
.linkBlackSmall:active
{
	color: #9b0521;
}
.linkBlackSmallBold
{
	font-weight:bold;	
	font-size: small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.linkBlackSmallBold:hover
{
	cursor: pointer;
	cursor: hand;
	color: #9b0521;
}
.linkBlackSmallBold:active
{
	color: #9b0521;
}
.linkSearchHitXSmall
{
	font-size: x-small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.linkSearchHitXSmall:hover
{
	cursor: pointer;
	cursor: hand;
	color: #cc9900;
}
.linkSearchHitXSmall:active
{
	color: #cc9900;
}
a.disabled
{
    pointer-events: none !important;
    cursor: default;
    color: #b3b3b3
}
a.special
{
	right: 25px;
	left: 386px;
	padding: 0px;
	width: 150px;
	position: relative;
	top: -4px;
	text-align: right;
}
a.what-is-orcid 
{	
	width: 100px;
	left: 185px;
	display: block;
	position: relative;
	text-align: left;
}
a.what-is-orcid-login 
{	
	width: 100px;
	display: block;
	position: relative;
	text-align: left;
}
.orcidRegButton 
{
    background-color:#cecece;
    background-image: linear-gradient(#ffffff 50%, #cecece 100%);
    background-image:-webkit-linear-gradient(top, #ffffff 50%, #cecece 100%);
    background-image:-moz-linear-gradient(top, #ffffff 50%, #cecece 100%);
    border:1px solid #cecece;
    padding:1px 5px;
    border-radius:0.35em;
    box-shadow: 1px 1px 2px #555555;
    cursor: pointer;
    color:#003365;
    line-height: 24px;
    vertical-align: middle;
}
.orcidRegButton:hover
{
    background-image: linear-gradient(#dedede 50%, #adadad 100%);
    border:1px solid #777777;
    cursor: pointer;
}
.orcidRegButton:active
{
    background-color: #e4e4e4;
    background-image: none;
    box-shadow: inset 1px 1px 2px #555555;
    outline: 0 !important;
}
.orcidRegButtonText 
{
    text-align:center;
    vertical-align: super;
    padding-left:5px;
    font-size: 85%;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.links
{
	margin-top: 10px;
	font-size: 100%;  /* CHECK THIS */
	left: 60px;
	margin-bottom: 10px;
	width: 400px;
	position: relative;
	text-align: center;
}
.links2
{
	margin-top: 10px;
	margin-bottom: 10px;
	text-align: center;
}
.links3
{
	margin-top: 10px;
	margin-bottom: 10px;
	margin-left: auto;
	margin-right: auto;
	font-size: 100%;
	text-align: center;
}

/*****************************************************************
	Styles for buttons
*****************************************************************/
div.buttons 
{
    text-align:center
}
div.buttons input 
{
    margin-left:1.5em;
}
div.buttons>:first-child 
{
     margin-left:0em
}
div.buttonAlign
{
	text-align: center;
	margin-top: 15px;
	margin-bottom: 5px;
}
div.buttonAlign input[type="button"], div.buttonAlign input[type="submit"], div.buttonAlign button[type="button"]
{
    margin: 0.5em;
}
div.buttonAlignEM
{
	text-align: center;
	padding-left: 1em;
	padding-right:1em;
	padding-top: 1em;
	padding-bottom: 1em;
}

/*****************************************************************
	Style for page titles
*****************************************************************/
h3
{
	padding-bottom: 4px;
	font-size: 150%;
}
.faqtitle
{
	font-weight: bold;
	font-size: 160%;
}
.bigtitle
{
	font-weight: bold;
	font-size: 140%;
}

/*****************************************************************
	Styles for various general text elements
*****************************************************************/
div.warningbox_enclosed
{
	margin-top: 0px;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 10px;
	border-top: #cc0000 1px solid;
	border-bottom: #cc0000 1px solid;
	border-right: #cc0000 1px solid;
	border-left: #cc0000 1px solid;
	background-color: #f4f4f4;
}
div.warningbox_enclosed_bold
{
	margin-top: 0px;
	padding-right: 10px;
	padding-left: 10px;
	padding-bottom: 10px;
	padding-top: 10px;
	border-top: #cc0000 1px solid;
	border-bottom: #cc0000 1px solid;
	border-right: #cc0000 1px solid;
	border-left: #cc0000 1px solid;
	background-color: #f4f4f4;
}
div.warningbox_enclosed_medium
{
	margin-top: 1em;
	padding-right: 1.5em;
	padding-left: 1.5em;
	padding-bottom: 1.5em;
	padding-top: 1.5em;
	margin-left: 2em;
	margin-right: 2em;
	border-top: #cc0000 double;
	border-bottom: #cc0000 double;
	border-right: #cc0000 double;
	border-left: #cc0000 double;
	background-color: #f4f4f4;
}
div.warningbox
{
	padding: 10px;
	margin-top: 20px;
	margin-bottom: 10px;
	margin-left: 10px;
	margin-right: 10px;
	border-top: #cc0000 1px solid;
	border-bottom: #cc0000 1px solid;
	left: 10px;
	width: 95%;
	background-color: #f4f4f4;
}
div.box
{
	padding: 10px;
	margin-top: 20px;
	margin-bottom: 10px;
	margin-left: 10px;
	margin-right: 10px;
	border-top: #002244 1px solid;
	border-bottom: #002244 1px solid;
	left: 10px;
	width: 95%;
	background-color: #f4f4f4;
}
div.boxClosed
{
	padding: 10px;
	margin-top: 20px;
	margin-bottom: 10px;
	margin-left: 10px;
	margin-right: 10px;
	border: #002244 1px solid;
	left: 10px;
	width: 95%;
	background-color: #f4f4f4;
}
.warning
{
	color: #cc0000;
}
.disabled-text
{
    color: #999999;
}
.errorText, .error-text
{
    font-weight: bold;
    color: #ff0000;
}
.redtext
{
	color: #cc0000;
}
.redtextBold
{
	color: #cc0000;
	font-weight: bold;
}
.greentext
{  
	color: #008800;
}
.orangetext
{
	color: #EE6600
}
.orangetext2
{
	color: #FF6600
}
.blacktext
{
	color: #000000;
}
.important2
{
	font-size: 105%;
}
.important2 a:link
{
	color: #003365;
	text-decoration: underline;
}
.important2 a:hover
{
	color: #6699cc;
}
.important2 a:active
{
	color: #ff3333;
}
.important3
{
	font-weight: bold;
}
.important3 a:link
{
	color: #003365;
	text-decoration: none;
}
.important3 a:hover
{
	color: #6699cc;
}
.important3 a:active
{
	color: #ff3333;
}
.important3 a:visited
{
	color: #003365;
	text-decoration: none;
}
.texthighlight
{
	font-weight: bold;
}
.texthighlightCentered
{
	font-weight: bold;
	text-align: center;
}
.important
{
	color: #003365;
	font-weight: bold;
}
p.content
{
	text-align:justify;
}
.centered
{
	text-align: center;
}
div.centered-div
{
    margin-left: auto;
    margin-right: auto;
}
.rightAligned
{
	text-align: right;
}
.leftAligned
{
	text-align: left;
}
.topAligned
{
    vertical-align: top;
}
.right
{
    float: right;
}
.alertText
{
	color: #cc0000;
}
.alertTextCentered
{
	color: #cc0000;
	text-align: center;
}
.italic
{
	font-style: italic;
}
.italic90
{
	font-style: italic;
    font-size: 90%;
}
.hiddenText
{
	color: #808080;
	font-style: italic;
}
.hide 
{
    display: none;
}

/*****************************************************************
	Styles for various general form elements
*****************************************************************/

/*
.faded_for_hidden_item
{
	filter:alpha(opacity=40);
	-moz-opacity:.40;
	opacity:.40;
}
*/
.faded_for_hidden_item
{
	color:#AAAAAA;
}
input.txt
{
	font-size: 100%;
	height: 18px;
	vertical-align: middle;
	line-height: 18px;
	width: 200px;
	color: #000000;
	margin-right: 5px;
	margin-top: 2px;
	margin-bottom: 2px;
	background-color: #ffffff;
	padding-top: 0px;
	padding-bottom: 0px;
}
input[type="button"].txt_disabled, input[type="button"].req_disabled
{
	background-color: #e0e0e0;
    background-image: none;
    box-shadow: none;
    font-size: 100%;
	height: 18px;
	vertical-align: middle;
	line-height: 18px;
    cursor: pointer;
	width: 206px;
	color: #000000;
	margin-right: 5px;
	margin-top: 2px;
	margin-bottom: 2px;
	padding-top: 0px;
	padding-bottom: 0px;
}
input.txt.disabled
{
	background-color: #e0e0e0;
}
.txtarea:disabled
{
	background-color: #e0e0e0;
}
.txtarea
{
	vertical-align: middle;
	color: #003366;
	background-color: #ffffff;
	padding: 0px;
	width: 100%;
}
/* class txtarea that ranges from 5 to 20 lines */
.txtarea_block
{
	border: #999999 1px solid;
	vertical-align: middle;
	color: #003366;
	background-color: #ffffff;
	padding-top: 0px;
	padding-bottom: 0px;
	width: 100%;
    min-height: 75px;
    // max-height is applied through javacript on EMDetails.aspx // TT 25624
}
/* class grey txtarea that ranges from 1 to 20 lines */
.txtarea_block_grey
{
	border: #f4f4f4 1px solid;
	vertical-align: middle;
	color: #003366;
	background-color: #f4f4f4;
	padding-top: 0px;
	padding-bottom: 0px;
	width: 100%;
    min-height: 15px;
    max-height: 280px;
}
.radio
{
	margin-top: 2px;
	margin-bottom: 2px;
	margin-left: 2px;
	vertical-align: middle;
}
.listmenu
{
	border: #003365 1px solid;
	vertical-align: middle;
	color: #003366;
	background-color: #ffffff;
	padding-top: 0px;
	padding-bottom: 0px;
	width: 225px;
}
.checkbox 
{ 
    vertical-align: middle; 
}
.disabledDisplay,
input.txt.disabledDisplay,
input.req.disabledDisplay,
input.radio.disabledDisplay,
input.button.disabledDisplay,
select.dropdown.disabledDisplay
{
    color: #000000;
    background-color: #e0e0e0;
    cursor: text;
}

/*****************************************************************
	Styles for general table data, including pagination & sorting
*****************************************************************/
td.leftWidth
{
	width: 22%; 
	background-color: #EEEEEE; 
}
td.leftColumn
{
	background-color: #EEEEEE; 
}
table.minitable 
{
	border:0;
	border-collapse: collapse;
	vertical-align: middle;
}
table.minitable td
{ 
	border:0;
	border-collapse: collapse;
}
table.minitable tr
{ 
	border:0;
	border-collapse: collapse;
}
th
{
	vertical-align: bottom;
}
.results
{
	padding: 2px;
	margin-top:5px;
	margin-bottom: 15px;
	margin-left: 0px;
	width: 796px;
	height: inherit;
	color: #003365;
	clear: both;
}
.reportResultsTable
{
	border: 0px;
	border-collapse:collapse;
	width: 100%;
	line-height: 14px;
	background-color: #FFFFFF
}
.reportResultsTable th
{
 	background-color:#ffffff;
	font-weight: bold;
	text-align: center;
	vertical-align:bottom;
	border: 0px;
}
.reportResultsTable td
{
	padding: 4px;
	border: 0px;
    vertical-align:top;
}
.reportResultsTable tr
{
	background-color: #EEEEEE;
}
.reportResultsTable tr.altrow
{
	background-color: #ffffff;
}
.datatable3
{
	border-collapse: collapse;
	width: 100%;
	line-height: 14px;
	background-color: #FFFFFF
}
.datatable3 td
{
	padding: 3px;
	border: 1px solid #ffffff;
}
.datatable3 th
{
	padding: 3px;
}
.datatable3 tr
{
	background-color: #f4f4f4;
}
.datatable3 tr.altrow
{
	background-color: #ffffff;
}
.datatable3 th
{
	border: 1px solid #ffffff;
	border-collapse: collapse;
	background-color: #3a5568;
	font-weight: bold;
	text-align: left;
	padding-left: 4px;
	color: #ffffff;
}
.datatable
{
	border: 1px solid #cccccc;
	border-collapse: collapse;
	width: 100%;
	line-height: 14px;
}
.datatable td
{
	padding: 3px;
	border: 1px solid #cccccc;
	vertical-align: top;
}
.datatable th
{
	padding: 3px;
}
.datatable tr
{
	background-color: #ffffff;
}
.autoRow tr:nth-child(even)
{
	background-color: #ffffff;
}
.autoRow tr:nth-child(odd)
{
	background-color: #e6e6e6;
}
.datatable th
{
	border-left: 1px #cccccc solid;
	background-color: #3a5568;
	font-weight: bold;
	text-align: left;
	padding-left: 4px;
	color: #ffffff;
}
.datatable caption
{
	color: #33517A;
	text-align: left;
	padding-top: 3px;
	padding-bottom: 8px;
}
.datatable tr.altrow
{
	background-color: #e6e6e6;
}
.datatable div
{
    padding: 0px;
}
.datatable2
{
	border: 1px solid #cccccc;
	border-collapse: collapse;
	width: 100%;
	line-height: 14px;
}
.datatable2 tr
{
	background-color: #ffffff;
	border: 1px solid #cccccc;
}
.datatable2 td
{
	border: 1px solid #cccccc;
	padding: 1px;
	vertical-align: top;
}
.datatable2 td.label
{
	border: 1px solid #cccccc;
	font-weight: bold;
	padding: 1px;
	vertical-align: top;
}
.datatable2 th
{
	border: 1px solid #cccccc;
	border-collapse: collapse;
	background-color: #3a5568;
	font-weight: bold;
	text-align: left;
	padding-left: 4px;
	color: #ffffff;
}
.datatable2 th.centered
{
	border: 1px solid #cccccc;
	border-collapse: collapse;
	background-color: #3a5568;
	font-weight: bold;
	text-align: center;
	padding-left: 4px;
	color: #ffffff;
}
.datatable2 caption
{
	color: #33517A;
	text-align: left;
	padding-top: 3px;
	padding-bottom: 8px;
}
.datatable2 tr.altrow
{
	background-color: #e6e6e6;
}
.datatable2 .metaHead
{
    color: black;
    font-weight: normal;
    font-style: italic;
    background-color: inherit;
    padding-left: 0;
}
.searchCriteriaNET
{
	border-collapse: collapse;
	line-height: 14px;
	text-align:center;
	margin-left: auto; 
	margin-right: auto; 
    border: solid 1px #cccccc;
}
.datatable4
{
	width: 100%;
	line-height: 14px;
}
.datatable4 tr.altrow
{
	background-color: #e6e6e6;
}
.searchCriteriaNET td
{
	padding: 3px;
}
.searchCriteriaNET th
{
	padding: 3px;
}


.searchCriteriaNET tr
{
    background-color: #ffffff;
}
.searchCriteriaNET th
{
	border-left: 1px #cccccc solid;
	background-color: #3a5568;
	font-weight: bold;
	text-align: left;
	padding-left: 4px;
	color: #ffffff;
}
.searchCriteriaNET caption
{
	color: #33517A;
	text-align: left;
	padding-top: 3px;
	padding-bottom: 8px;
}
.searchCriteriaNET tr.altrow
{
	background-color: #e6e6e6;
}
.searchCriteriaSmallColumn 
{
    width: 1%;
}
.searchCriteriaNegation 
{
    width: 10%;
}
.searchCriteriaSelector 
{
    width: 20%;
}
.searchCriteriaValue 
{
    width: 20%;
}
.searchCriteriaValue 
{
    width: 25%;
}
.paginationTable
{
	border: 0px;
	border-collapse: collapse;
	width: 100%;
	line-height: 14px;
}
.paginationTable tr
{
	background-color: #f4f4f4;
}
.paginationTable td
{
	border: 0px;
	padding: 1px;
}
td.pagination
{
	text-align: center;
}

/*****************************************************************
	Styles for custom flags
*****************************************************************/
span.flag
{
	font-family: "Zapf Dingbat", "Arial Unicode MS", "Code2000", "Lucida Sans Unicode";
    font-size:175%;
}
span.disabledflag
{
	font-family: "Zapf Dingbat", "Arial Unicode MS", "Code2000", "Lucida Sans Unicode";
	font-size: 175%;
	color: #cccccc;
	line-height: 110%;
}
table.flagmatrix
{
	text-align: center;
	vertical-align: middle;			
	background: #ffffff;
	border: 0px #ffffff hidden;
}
table.flagmatrix td.disabledflag
{
	border: 2px #ffffff hidden;
	background: #dddddd;
}
table.flagmatrix td.selectedflag
{
	border: 2px #000000 dashed;
}
table.flagmatrix td.unselectedflag
{
	border: 2px #ffffff solid;
}	
div.collapsedMenuFlags
{
	white-space: normal;
}
div.flaggroup
{
	width: 75px; 
	white-space: normal;
}
div.actionMenuLinksWithFlags
{
	width: 200px; 
	white-space: normal;
	vertical-align: middle;
}

/*****************************************************************
	Miscellaneous other styles
*****************************************************************/
br
{
	clear: left;
}
img
{
	border: 0px;
}
img.calendar
{
	vertical-align: middle;
}
img.smallicon
{
    width:14px;
    height:14px;
    border:0px;
}
p.acrobatLink
{
	font-weight: bold;
	text-align: center;
}
th.shaded
{
	border-left: 1px #cccccc solid;
	background-color: #cccccc;
	font-weight: bold;
	text-align: left;
	padding-left: 3px;
	color: #ffffff;
}
td.shaded 
{
    background-color:#EEEEEE;
}
td.shadedNoWrap
{
    background-color:#dddddd;
    white-space:nowrap
}
td.shadedBold
{
    background-color:#EEEEEE;
    font-weight:bold;
}
span.proxy 
{ 
    background: #fcee00; 
}
.proxy 
{ 
    background: #fcee00; 
}
.nowrap 
{ 
    white-space:nowrap 
}
.lettercontent 
{
	margin-top:5px;
	border-top:1px solid #000033;
	border-bottom:1px solid #000033;
	padding:5px; 
	background-color:#FFFFFF;
    word-break: break-word;
}
/* Links on the top and the bottom of the checkbox column for selecting/deselecting all */
a.checkAll 
{
    white-space: nowrap;
    display: inline;
}

/*Spec 12.037 20150514 GBS*/
span.label_prxyselfreg
{
	padding-right: 2em;
	display: block;
	float: left;
	width: 165px;
	color: #003366;
	text-align: right;
    font-weight:bold;
	line-height: 20px;
	vertical-align: middle;
	padding-top: 0px;
	padding-bottom: 0px;	
}

/**********************************************************************************************************
********
********  FEATURE-SPECIFIC STYLES
********
***********************************************************************************************************/

/*****************************************************************
	Styles for registration, UMI, and UOPI pages. Form fields & labels
*****************************************************************/
span.labeltagshort
{
	float: left;
	text-align: left;
	padding-right: 2em;
	width: 100px;
	color:#003366;
	line-height: 20px;
	vertical-align: middle;
	padding-top: 0px;
	padding-bottom: 0px;
}
span.label
{
	display: block;
	float: left;
	text-align: right;
	padding-right: 2em;
	width: 165px;
	color:#003366;
	line-height: 20px;
	vertical-align: middle;
	padding-top: 0px;
	padding-bottom: 0px;	
}
span.labelLeft
{
	display: block;
	float: left;
	text-align: left;
	padding-right: 2em;
	width: 160px;
	color:#003366;
	line-height: 20px;
	vertical-align: middle;
	padding-top: 0px;
	padding-bottom: 0px;	
}
span.label_req
{
	padding-right: 2em;
	display: block;
	float: left;
	width: 165px;
	color: #cc0000;
	text-align: right;
	line-height: 20px;
	vertical-align: middle;
	padding-top: 0px;
	padding-bottom: 0px;	
}
.helpersreq
{
	font-size: 90%;
	color: #cc0000;
	text-align: left;
}
.helpers
{
	font-size: 90%;
	color: #003366;
	text-align: left;
}
.hint
{
	border-right: #333333 1px;
	border-top: #333333 1px;
	border-left: #333333 1px;
	border-bottom: #333333 1px;
	padding-left: 2px;
	padding-bottom: 2px;
	padding-top: 2px;
	padding-right: 2px;
	font-size: 90%;
	width: 320px;
	left: 185px;
	color: #003366;
	display: block;
	position: relative;
	background-color: #ffffff;
	text-align: left;
}
.hint_sample
{
	border-right: #333333 1px;
	border-top: #333333 1px;
	border-left: #333333 1px;
	border-bottom: #333333 1px;
	padding-left: 2px;
	padding-bottom: 2px;
	padding-top: 2px;
	padding-right: 2px;
	font-size: 90%;
	width: 320px;
	left: 185px;
	color: #003366;
	display: block;
	position: relative;
	text-align: left;
}
input.req
{
	font-size: 100%;
	height: 18px;
	vertical-align: middle;
	line-height: 18px;
	width: 200px;
	color: #000000;
	margin-right: 5px;
	margin-top: 2px;
	margin-bottom: 2px;	
	padding-top: 0px;
	padding-bottom: 0px;
}
.radioreq
{
	margin-top: 2px;
	margin-bottom: 2px;	
	margin-left: 2px;
	vertical-align: middle;
}

/*****************************************************************
	Styles for manuscript submission step menu
*****************************************************************/
ul.stepitems
{
	width: 184px;
	padding-left: 0px;
	padding-right: 10px;
	margin-left: 0px;
	list-style-type: none;
	display: inline;
	font-weight: bold;
	font-size: 90%;
	text-decoration: none;
	list-style-position: outside;
	margin-top: 5px;
	margin-bottom: 5px;
}

ul.stepitems li
{
	padding: 0px 0px 0px 0px;
	border: 1px solid #003366;
	line-height: 170%;
	list-style-image: url(img/blank-bullet.gif);
	list-style-position: outside;
	margin-top: 5px;
	margin-bottom: 5px;
	margin-left: 15px;
	background-color: #F4F4F4;
	border:none;
	padding:0.25em 0.5em;;
	width:180px;
}
ul.stepitems li.check
{
	border: 0px solid #003366;
	list-style-image: url(./img/checkmark.gif);
	color: #cccccc;
}
ul.stepitems li.arrow
{
	list-style-image: url(./img/arrow.gif);
	color: #002244;
	background-color: #cccccc;
}
/* Added Spec 13.0-36 ALAW 20151203*/
ul.stepitems li.warning
{
	list-style-image: url(./img/doc_status_warning_sm.png);
	color: #cccccc;
}
ul.stepitems li a
{
	display: block;
	color: #002244;
	text-decoration: none;
}
ul.stepitems li a:visited
{
	display: block;
	color: #002244;
	text-decoration: none;
}
ul.stepitems li:hover
{
	background-color:#cccccc;
	color: #002244;
	text-decoration: none;
}
ul.stepitems li a:active
{
	display: block;
	background-color:#cccccc;
	color: #002233;
	text-decoration: none;
}

/*****************************************************************
	Styles for manuscript submission step menu
	for "Short" submission interface. (Spec 7.0-38 DUS 20090220)
*****************************************************************/
div.stepitems
{
	width: 184px;
	padding-left: 0px;
	padding-right: 10px;
	margin-left: 0px;
	list-style-type: none;
	font-weight: bold;
	font-size: 90%;
	text-decoration: none;
	list-style-position: outside;
	margin-top: 5px;
	margin-bottom: 5px;
}

div.stepitems div
{
	border: 1px solid #003366;
	line-height: 170%;
	list-style-image: url(img/blank-bullet.gif);
	list-style-position: outside;
	margin-top: 5px;
	margin-bottom: 5px;
	margin-left: 15px;
	background-color: #F4F4F4;
	border:none;
	width:190px;
}
div.stepitems div.currentStep
{
	list-style-image: url(./img/arrow.gif);
	background-color: #cccccc;
 	width:180px;
    padding:0.25em 0.5em;
}
div.stepitems div.doneStep
{
	list-style-image: url(./img/checkmark.gif);
	width:190px;
}
div.stepitems a
{
	display: block;
	color: #002244;
	text-decoration: none;
	width:180px;
	padding:0.25em 0.5em;;
}
div.stepitems a:visited
{
	display: block;
	width:180px;
	color: #002244;	
 }
div.stepitems a:hover
{
	background-color:#cccccc;
	color: #002244;
	width:180px;
}
div.stepitems a:active
{
 	background-color:#cccccc;
	color: #002233;
	width:180px;
}

/*****************************************************************
	Styles for Login pages
*****************************************************************/
.loginError
{
	color: #cc0000;
	font-weight: bold;
}
div.copyright
{
	text-align: center;
	font-size: 85%;
}
div.copyright span
{
    padding: 0 1em 0 1em;
}
div.journalThumbnail
{
	text-align: center;
}

/* ****************************************************************
	Styles for reviewer & editor recommendation pages
*****************************************************************/
input.rating
{
	border: #003365 1px solid;
	margin-top: -10px;
	height: 18px;
	vertical-align: middle;
	line-height: 18px;
	width: 30px;
	color: #003366;
	margin-right: 5px;
	background-color: #ffffff;
	padding-top: 0px;
	padding-bottom: 0px;
}

/*****************************************************************
	Styles for classifications
*****************************************************************/
td.classificationscell
{ 
	padding: 2px; 
	border-bottom: 1px solid #000000; 
	background-color:#ffffff;
}
td.classificationsshadedcell
{ 
	padding: 2px; 
	border-bottom: 1px solid #000000; 
	background-color: #cccccc;
}
th.classifications
{ 
	padding: 0px; 
	border-bottom: 1px solid #000000; 
	background-color:#ffffff;
}
table.classifications
{
	width: 100%;
	border-collapse: collapse;
	border: 1px solid #999999;
	background-color: #000000;
	list-style-position: outside;
}

/*****************************************************************
	Style for ident.asp
*****************************************************************/
.navbarText
{
	color: #ffffff; 
	padding: 0px 0px 0px 0px;
	margin: 0px 0px 0px 0px;
}

/*****************************************************************
	Styles for Main Menus
*****************************************************************/
.main_menu_item
{
	line-height: 180%; 
	padding-left: 100px;
}
.main_menu_item_heading
{
	margin-left:-30px;
	margin-top:-20px;
}
.main_menu_item_heading_sub
{
    padding-left:2em;
    display:block;
}
#header 
{
	float:left;
	width:100%;
	background: url("../pone/img/bg-blackbar.gif") repeat-x bottom;
	font-size:93%;
	line-height:normal;
}
#header ul 
{
	margin-top:0;
	padding:0 10px;
	list-style:none;
}	
#header li
{
	float:left;
	background:#336699 url("../pone/img/norm_left2.gif") no-repeat left top;
	margin:0;
	padding:0 0 0 9px;
}
#header a 
{
	float:left;
	display:block;
	background:url("../pone/img/norm_right2.gif") no-repeat right top;
	padding:5px 15px 4px 6px;
	text-decoration:none;
	font-weight:bold;
	color:#FFFFFF;
}
#header a:hover 
{
	color:#FFFFFF;
	text-decoration: underline;
}
#header #current 
{
	background:#FFFFFF url("../pone/img/norm_left2.gif") no-repeat left top;
}
#header #current a 
{
	background:url("../pone/img/norm_right2.gif") no-repeat right top;
	color:#336699;
	padding-bottom:5px;
}
#header #current a:hover
{
	text-decoration: none;
}
div.menuSection
{
	text-align:center;
	border:1px solid #CCCCCC;
	/*font-size:10px;*/
}

/* These are the style specifications for the menu layers. */
div.menu  
{ 
	position:absolute; 
	visibility:hidden; 
	left:0; top:0;
	z-index:500;
	margin: 0px;
	padding:5px;
	border: 3px solid;	
	border-color: #CCC #666 #666 #CCC;
	background-color:#DDDDDD; 
	text-align:left;
	white-space:normal; /**/
}
div.menu a 
{ 
	cursor: pointer;
	padding:.25em .4em; 
	margin:0; 
	background-color:transparent; 
	display:block; 
	position:relative; 
	text-decoration:none; 
	line-height:11px;
}
div.menu a:link 
{ 
    text-decoration:none; 
}
div.menu a:visited 
{
    text-decoration:none; 
}
div.menu a:hover 
{ 
    color:#fff; 
    background-color:#3a5568; 
    text-decoration:none; 
}
/* Special action link styles for action links combined with flags */
div.menu a.linkWithFlags
{ 
	display: inline;
}
.false_link
{
	color: #0000ff;
	text-decoration: underline;
	cursor: pointer;
}
div.menu-expand a 
{ 
	cursor: pointer;
	padding:.25em .4em; 
	margin:0; 
	background-color:transparent; 
	display:block; 
	position:relative; 
	line-height:11px;
}

/********************************************
    Miscellaneous Font Styles
********************************************/
.normalWeight
{
    font-weight:normal;
}
.lblErrorMsgSmall
{
	font-size: x-small;
	color: #ff0000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblErrorMsgSmallBold
{
	font-size: x-small;
	color: #ff0000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.lblErrorMsgSmallBoldItalic
{
	font-size: x-small;
	color: #ff0000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
	font-style: italic;
}
.lblErrorMsgMedium			
{
	font-size: small;
	color: red;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblErrorMsgLarge
{
	font-size: medium;  
	color: red;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblInfoMsgSmall
{
	font-size: x-small;
	color: #003399;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblInfoMsgSmallBold
{
	font-size: x-small;
	color: #003399;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.lblInfoMsgMedium
{
	font-size: small;
	color: #003399;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblInfoMsgMediumBold
{
	font-size: small;
	color: #003399;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.lblInfoMsgLarge
{
	font-size: medium;
	color: #003399;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblSearchHitXSmall
{
	font-size: x-small;
	color: #000000;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.lblMissingItemSmall
{
	font-size: x-small;
	color:#AAAAAA;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblMissingItemSmallBold
{
	font-size: x-small;
	color:#AAAAAA;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.lblInactiveMsgSmall
{
	font-size: x-small;
	color: #808080;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
}
.lblInactiveMsgSmallBold
{
	font-size: x-small;
	color: #808080;
	font-family: verdana, geneva, arial, helvetica, sans-serif;
	font-weight: bold;
}
.white_links a:link
{
	color: #FFFFFF;
}
.white_links a:visited
{
	color: #FFFFFF;
}
.white_links a:hover
{
	color: #FFFFFF;
}
.disabledText
{
	color: #999999; 
	font-style: italic;
}
.bold
{
	font-weight: bold;
}

.radioButtonControl td
{
    width:36px;
    text-align:center;
}
#mbox
{
	background-color:#eee;
	padding:8px; 
	border:2px outset #666;
}
#mbm
{
	font-family:sans-serif;
	font-weight:bold;
	float:right;
	padding-bottom:5px;
}
#ol
{
	background-image: url(overlay.png);
}
.dialog 
{
	display:none
}
.motionDialog
{
    border-style:solid;
    border-width:thin;
    position:absolute;
    background-color:#EBEBEB;
    cursor:default;
    left:0px;
    top:0px;
    display:none;
}
.motionPopUp
{
    border-style:solid;
    border-width:thin;
    position:absolute;
    background-color:white;
    cursor:default;
    left:0px;
    padding:10px;
    display:none;
}
.validationSummary ul
{
    display:none
}
.validationSummary 
{
    height:100%;
    font-size:10pt;
    font-family:arial,helvetica,sans-serif;
}
.textMenu
{
	border: 0;
	padding: 0;
	margin: 0;
	width: 100%;
	font-family: Trebuchet MS, verdana, geneva, arial, helvetica, sans-serif;
	font-size: 11px;
	line-height: 13px;
	height: 85px;
	background: url(/pone/img/maintile.jpg);
}
.textMenu > .menuImage
{
	height: 53px;
	overflow: hidden;
}
.textMenu > ul
{
	margin: 0;
	padding: 0;
	white-space: nowrap;
}
.textMenu > ul > li
{
	display: inline;
	padding-left: 11px;
	padding-right: 4px;
	text-transform: uppercase;
	background: url(/pone/img/menuBullet.gif) no-repeat left center;
}
.textMenu > ul > li:first-child
{
	background: none;
	padding-left: 3px;
	padding-right: 4px;
}
.textMenu > ul > li > a
{
	text-decoration: none;
	white-space: nowrap;
}
.textMenu > ul > li > a:hover
{
	text-decoration: underline;
}
.textMenu > ul > li > a:active
{
	color: #CCCCCC;
	text-decoration: underline;
}
div.SpecialCharacters
{	
    Z-index:1;	
    position:fixed;	
    max-width:55.25em;	
    right:1.5%;	top:1.5%;	
    background-color:#DAE3E7;	
    opacity:1.0;	
    filter:Alpha(opacity=100);	
    border:#AAAAAA 1px solid;	
    border-radius:0.5em;	
    padding:2%;	
     
    cursor: move
}
div.SpecialCharacters a
{	
    text-decoration:none;
	color:#003365;	
    display:inline-block;	
    padding:0.25em;
	font-weight:bold;
}
.numbers
{
	font-size: 140%;
 	font-family:Cambria,"Times New Roman",Serif;
	font-weight:bold;
}
div.MergeCodes
{
	Z-index:1;
    position:fixed;
	opacity:1.0;
    filter:Alpha(opacity=100);
	border:#AAAAAA 1px solid;
	border-radius:0.5em;
	padding:2%;
	font-size:110%;
    cursor: move
}
tr.gridRowStyle
{
    height: 2.5em;
}
tr.gridHeaderStyle
{
    height: 4em;
}
span.incompleteIndicator
{
	color: blue;
}

div.score
{
    padding-left:0.5em; 
}

.preLine > span
{
    white-space:pre-wrap;
}
.preLine p
{
    white-space:pre-wrap;
}

/*****************************************************************
	Styles for expandable text area
******************************************************************/

.expand, 
.expandingArea
{
    position: relative;
    border: 1px solid #ccc;
    font: 100% verdana, geneva, arial, helvetica, sans-serif;
    width: 99%;
}
       
.unavailable 
{
    color: gray;
    font-style: italic;
}
.clear 
{
    clear: both;
}

/*****************************************************************
	Styles for ActionLock
*****************************************************************/
.lockMessageDiv 
{
    height: 40px;
    width: 60%;
    margin: auto;
    font-weight: bold;
    text-align:center;
    line-height: 40px;
    background-color: #FFFFCC;
    border: 1px solid #000000;
    margin-bottom: 1em;
}
.availableMessageDiv 
{
    height: 40px;
    width: 60%;
    margin: auto;
    font-weight: bold;
    text-align:center;
    line-height: 40px;
    background-color: #CCFF99;
    border: 1px solid #000000;
    margin-bottom: 1em;
}
.blockDiv 
{
    margin: 15px 5px 5px 5px;
    padding:  5px 5px 5px 5px;
}

/*****************************************************************
	Styles for calendar
*****************************************************************/
.cal_Theme .ajax__calendar_container   
{
    background-color: #4682B4;
    border:solid 1px #4682B4;
    padding: 1px;
}
.cal_Theme .ajax__calendar_header  
{
    background-color: #4682B4; 
}
.cal_Theme .ajax__calendar_title,
.cal_Theme .ajax__calendar_next,
.cal_Theme .ajax__calendar_prev    
{
    color: #ffffff;
    font-weight: normal;
    font-size: small;
    font-family: tahoma, verdana;
    padding-top: 3px;
}
.cal_Theme .ajax__calendar_prev 
{
    background-image: url(./img/prev.gif);
    width: 16px;
    height: 16px;
}
.cal_Theme .ajax__calendar_next 
{
    background-image: url(./img/next.gif);
    width: 16px;
    height: 16px;
}
.cal_Theme .ajax__calendar_body    
{
    text-align: center;
    padding-bottom: 2px;
    margin: 0px;
}
.cal_Theme .ajax__calendar_days,
.cal_Theme .ajax__calendar_days td  
{
    background-color: #87CEFA;
    border: solid 1px #4682B4;
    text-align:center;
    padding: 0px;
    margin: 0px;
    border-left: 0px;
}
.cal_Theme .ajax__calendar_days tr,
.cal_Theme .ajax__calendar_days table,
.cal_Theme .ajax__calendar_days div 
{
    border: 0px;
}
.cal_Theme .ajax__calendar_days table,
.cal_Theme .ajax__calendar_months table,
.cal_Theme .ajax__calendar_years table  
{
    width: 100%;
    height: 100%;
    margin: 0px;
}
.cal_Theme .ajax__calendar_dayname
 {
    background-color: #87CEFA;
    color: #ffffff;
    text-align: center;
    height: 18px;
    width: 19px;
}
.cal_Theme .ajax__calendar_day
{
    background-color: #FFFFFF;
    text-align: right;
    text-decoration: underline;
    text-decoration-color: blue;
    height: 18px;
    width: 19px;
}
.cal_Theme .ajax__calendar_invalid .ajax__calendar_day  
{
    text-decoration: none;
}
.cal_Theme .ajax__calendar_invalid.ajax__calendar_hover .ajax__calendar_day  
{
    font-weight: normal;
    cursor: default;
}
.cal_Theme .ajax__calendar_hover .ajax__calendar_day,
.cal_Theme .ajax__calendar_hover .ajax__calendar_month,
.cal_Theme .ajax__calendar_hover .ajax__calendar_year  
{
    color: #004080; 
    font-weight: bold; 
    background-color: #ffffff;
}
.cal_Theme .ajax__calendar_active .ajax__calendar_day,
.cal_Theme .ajax__calendar_active .ajax__calendar_month,
.cal_Theme .ajax__calendar_active .ajax__calendar_year  
{
    color: #004080; 
    font-weight: bold; 
    background-color: #ffb6c1;
}
.cal_Theme .ajax__calendar_footer   
{
    display: none;
}
.cal_Theme .ajax__calendar_today  
{
    font-weight:normal;
    background-color: #ffb6c1;
}
.cal_Theme .ajax__calendar_other,
.cal_Theme .ajax__calendar_hover .ajax__calendar_today,
.cal_Theme .ajax__calendar_hover .ajax__calendar_title 
{
    color: #bbbbbb;
}

/*****************************************************************
	Styles for Auto Suggest control
*****************************************************************/
 /* http://docs.jquery.com/UI/Autocomplete#theming */
.ui-autocomplete 
{ 
    position:absolute; 
    cursor:default; 
    background:#CCCCC;
}   

/*****************************************************************
	Styles for admin pages
*****************************************************************/
th.pageText, td.pageText, p.pageText, div.pageText
{
	font-size: small;
}


span.vertical-validator[style*="inline"]
{
	display:block !Important;
}
span.horizontal-validator
{
	padding-right: 0.5em;
}
/* Overlay and Loader styles */
.overlay
{
    position: fixed;
    z-index: 99;
    top: 0px;
    left: 0px;
    width: 100%;
    height: 100%;
    background-color: #000;
	filter: Alpha(Opacity=0.2);
	opacity: 0.2;
	-moz-opacity: 0.2;
}
* html .overlay
{
    position: absolute;
    height: expression(document.body.scrollHeight > document.body.offsetHeight ? document.body.scrollHeight : document.body.offsetHeight + 'px');
    width: expression(document.body.scrollWidth > document.body.offsetWidth ? document.body.scrollWidth : document.body.offsetWidth + 'px');
}
.loader
{
    z-index: 100;
    position: fixed;
    width: 120px;
    margin-left: -60px;
    top: 50%;
    left: 50%; 
}
* html .loader
{
    position: absolute;
    margin-top: expression((document.body.scrollHeight / 4) + (0 - parseInt(this.offsetParent.clientHeight / 2) + (document.documentElement && document.documentElement.scrollTop || document.body.scrollTop)) + 'px');
}  
/*	Spec 11.0-04 */
.ui-autocomplete  
{
	position:absolute; 
	cursor:default; 
	background:#CCCCC;
} 
.ui-autocomplete-input
{
	height: 18px;
	font-size: 11.2px;
}
.custom-combobox 
{
	position: relative;
	display: inline-block;
	height: 18px;
    width: 207px;
	font-size: 11px;
	line-height: 13px;
    white-space: nowrap;
}
.custom-combobox a.custom-combobox-toggle
{
	position:relative;
	top: 0;
	bottom: 0;
    margin-top: 1px;
    margin-bottom: 1px;
	margin-left: -1px;
	padding: 0;
    /* support: IE7 */
	height: 18px; /* 1.8em; */
    width: 24px;
	font-size: 11px;
	line-height: 18px;
    border-width: 1px;
    border-color: #f0f0f0;
    border-style: inset;
    border-radius: 0px;
}
input.custom-combobox-input,  select.custom-combobox-input
{
	/* margin: 0; */
	padding: 0px 0px 0px 2px;  /* 0.3em; */
	width: 150px;
	height: 18px;
	/* font: MS Shell Dlg; */
	font-size: 11px;
	color: rgb(0, 0, 0);
	margin-right: -1px;
	margin-top: 1px;
	margin-bottom: 1px;
	background-color: rgb(255, 255, 255);
	vertical-align: middle;
	line-height: 18px;
    border-width: 1px;
    border-color: #f0f0f0;
    border-style: inset;
    border-radius: 0px;
}
.autocompleteinput
{
	font-size: 100%;
    border: 1px solid rgb(153, 153, 153);
    height: 18px;
    vertical-align: middle;
    line-height: 18px;
    width: 200px;
    color: rgb(0, 0, 0);
    margin-right: 5px;
    margin-top: 2px;
    margin-bottom: 2px;
    background-color: rgb(255, 255, 255);
    padding-top: 0px;
    padding-bottom: 0px;
}
/*  Spec #: 11.0-26 styles for aries.selectable, aries.sortable widgets */
.datatable tr.em-ui-selected, .datatable tr.em-ui-selected td
{
    background-color: #f4f4f4;
}
.em-ui-state-default>td:first-child 
{
    background-image: url('img/drag.png');  /* BUGFIX 27437 20140827 BBD */
    background-position: center left;       /* BUGFIX 27437 20140827 BBD */
    background-repeat: no-repeat;
    padding-left:20px;
}
.em-ui-sortable-placeholder
{
    height: 3px !important;
}        
.em-ui-sortable-placeholder td
{
    line-height: 3px;
    border-top: 1px solid gray;
    border-bottom: 1px solid gray;
    background-color: gray;
    *padding:0; /* IE 7 and below */
}
.em-ui-sortable-placeholder td hr
{
    border-top: 1px dashed white;
    border-bottom: none;
    margin:0;
    display:block;
    *margin: -7px 0 0 0;    /* IE 7 and below */
    *display: list-item;    /* IE 7 and below */
    padding-left: 5px;
    padding-right: 5px;
}
.loading
{
    background-image: url('img/AjaxLoaderBlue40.gif');
    background-position: center center;
    background-repeat: no-repeat;
    min-height: 100px;
}
.spinner-wrapper 
{
    background-color: white;
    padding: 1.0em;
    opacity: 0.9;
    display:table;
    border-radius:0.35em;
    -webkit-border-radius:0.35em;
    -moz-border-radius:0.35em;
    box-shadow: 1px 1px 2px #555555;
    margin: auto;
}
.spinner
{
    margin: auto;
    background-image: url('img/AjaxLoaderBlue40.gif');
    background-position: center;
    background-repeat: no-repeat;
    width: 40px;
    height: 40px;
}
.spinner-message 
{
    margin-top: 0.5em;
    color: #336699;
    font-weight: bold;
    text-align:center;
}
.viewContent 
{
    text-align: center;
}
.highlightBorder
{
    border-color: #FF9933 !important;
}

.status-message 
{
	margin: 0;
	padding: 8px;
	background: #f4f4f4;
	border-radius: 8px;
	margin:10px 0 10px 10px;
}

.AboutLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/about.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}

.TotalsLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/totals.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}
.DownloadLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/download.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}
.CancelLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;     
    background-image: url(img/cancel.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}
.SearchLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/search.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}
.SaveLink 
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/save.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
    color: transparent;
}
.CheckAllLink 
{
    height: 22px;
    width: 22px;
    background-image: url(img/checkall.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
}
.ClearAllLink 
{
    height: 20px;
    width: 20px;
    background-image: url(img/clearall.png);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
}
.ToggleLink, #ToggleLink
{
    height: 22px;
    width: 22px;
    margin: 0 0.75em 0 0;    
    background-image: url(img/toggle.png?v=1);
    background-repeat: no-repeat;
    display: inline-block;
    background-size: contain;
}
.linkStyle 
{
    cursor:pointer;
    color:blue;
}
.linkColor
{
    color:blue;
}
.ro-overlay-content 
{
    position: absolute;
    z-index:101;
    border: 1px solid #fff; 
    background-color: #f4f4f4;
    color: #333;
    top: 0px;
    left: 0px;
}
.ro-hidden 
{
    display: none;
}
.ro-empty 
{
}
.ro-container 
{
    height: 100%;
    width: 100%;
    display: table-cell;
    vertical-align: middle;
}
.reviewerHoverDetails 
{
    margin: auto 0;
    display: inline-block;
}
.ui-icon-pin-w, .ui-icon-pin-s 
{
    display: inline-block;
    position: absolute;
    right: 0;
}
.ui-icon-pin-w:hover, .ui-icon-pin-s:hover 
{
    background-color: #cccccc;
}
div.ro-content 
{
    height: 100%;
}
div.ro-loading 
{
    text-align: center;
    background-color: #f4f4f4;
    height: 100%;
    font: 0/0 a;
}
div.ro-loading:before 
{
    content: ' ';
    display: inline-block;
    vertical-align: middle;
    height: 100%;
}
div.ro-loading > img 
{
    display: inline-block;
    vertical-align: middle;
}
div.rs-less-link,
div.rs-more-link 
{
    display:table-caption;
    caption-side:bottom;
    padding-top:0.5em;
}
.ro-overlay-content .rs-less-link 
{
    display: none;
}
.popup-shadow 
{
    background-color:#ffffff;
    border-radius:0.5em;
    box-shadow: 0 0 10px #336699;
}
.ro-overlay 
{
}
.rs-header 
{
    display:block;
    white-space:nowrap;
}
.rs-complete 
{
	color: #008800;
}
.rs-partialReview  
{
    color:#009999;    
}
.rs-agreed  
{
    color:#EE6600;             
}
.rs-invited 
{
    color: #FF0000;
}
.rs-declined 
{
}
.rs-late 
{
    color:#FF0000;
}
div.rs-table, 
div.roi-table 
{
    display: table;
}
div.rs-headerRow
{
    display: table-row;
}
div.roi-table > .rs-headerRow:not(:first-child)
{
    height: 1.5em;
}
div.roi-table > .rs-headerRow:not(:first-child) .rs-countCell
{
    padding-top : 0.5em;
}
div.rs-detailCell
{
    display: table-cell;
    white-space:nowrap;
}
.roi-table >.roi-detailRow 
{
    display: table-row;
    font-size: smaller;
    font-style: italic; 
}
.roi-table > caption
{
    border-top:1em;
}
.rs-headerRow  >.rs-detailCell,
.roi-detailRow  >.roi-detailCell   
{
    display: table-cell;
} 
.rs-headerRow >.rs-countCell  
{
    display: table-cell; 
    padding-right: 0.5em;
}
.roi-detailRow >.roi-countCell,  
{
    display: table-cell; 
    padding-right: 0.5em;
}

.roi-detailRow >.roi-countTotalRequired,
.roi-detailRow >.roi-elapseDate,
.roi-detailRow >.roi-suggested,
.roi-detailRow >.roi-duedate,
.roi-grey 
{
    color: grey;
}
.roi-infoImage 
{
    vertical-align:middle; 
    margin-left:2px;
    margin-bottom:2px; 
    width:10px; 
    height:10px
}
.roi-overdue 
{
    color: #FF0000;       
}
.normalFont
{
    font-style: normal;
}
.roi-reasons-container 
{
    max-height: 300px;
    overflow: auto;
}
.ro-reviewer-reason-row 
{
    background-color: white;
}
.ro-reviewer-reason-row:not(:last-child) 
{
    padding-bottom: 1em;
}
.ro-reviewer-reason-name 
{
    padding-bottom: 0.5em;
    font-weight: bold;
}
.roi-view-reasons 
{
    font-size: 80%;
}
.DPNCDelete
{
    display: block;
    margin: 0 auto;
}
.DPNCDelete:hover
{
    cursor: pointer;
}
.DPNCSortButton:hover
{
    cursor: pointer;
}
.waitHover:hover
{
    cursor: wait;
}
.DPNCGridContainer
{
    clear:both;
}
.DPNCGridContainer.nonPrint
{
    height: 22em;
    overflow: scroll;
    resize: both;
}
.DPNCAlertText
{
    margin-top: 15px;
}
.DPNCRightSpacing
{
    margin-right: 1em;
}
.DPNCAlert.ui-dialog .ui-dialog-titlebar
{
    display:none;
}
.DPNCAlertContainer
{
    display:none;
}
.DPNCLabel
{
    color: blue;
    font-weight: bold;
    vertical-align: middle;
}
.DPNCLabel:hover
{
    cursor: pointer;
}

.error 
{
	font-size: small;
	color: #ff0000;
	font-weight: bold;
    padding-right: 2px;
}

.error_container 
{
	font-size: x-small;
	color: #ff0000;
	font-weight: bold;
    padding-top: 10px;
    padding-bottom: 10px;
}        

.hide 
{
    display:none;
}
.searchControl_ContainerTable
{
    margin-right: auto;
    margin-left: auto;
    border: 0;
}
.searchControl_ContainerTable table
{
    width: 781px;
}
.searchControl_CriteriaTable
{
    border-left: solid 1px #cccccc ;
    border-right: solid 1px #cccccc;
    border-top: solid 1px #cccccc;
    border-spacing: 0;
}
.searchControl_CriteriaTable td
{
    border-bottom: solid 1px #cccccc;
    padding-top: 4px;
    padding-bottom: 4px;
    padding-left: 2px;
    padding-right: 2px;
    background-color: #ffffff;
}
.searchControl_CriteriaTable th
{
    background-color: #3a5568;
    font-weight: bold;
    text-align: left;
    padding-left: 4px;
    color: #ffffff;
    vertical-align: bottom;
    border-right: solid 1px #ffffff;
}
.searchControl_ButtonRow
{
    width: 100%;
    border-left: solid 1px #cccccc;
    border-right: solid 1px #cccccc;
    border-bottom: solid 1px #cccccc;
    background-color: #ffffff;
}

.superScript
{
    vertical-align: super;
}
.strikeThrough
{
    text-decoration: line-through;
}

.cke_top 
{
    background-image: none !Important;
    background-color: #336699 !Important;
    background: -webkit-linear-gradient(90deg, #75a8da, #336699) !Important;
    background: -o-linear-gradient(0deg, #75a8da, #336699) !Important;
    background: -moz-linear-gradient(0deg, #75a8da, #336699) !Important;
    background: linear-gradient(0deg, #75a8da, #336699) !Important;
}

.cke_dialog_ui_vbox_child div 
{
    width:100% !Important;
}
.cke_pasteframe 
{
    width:100% !Important;
}
.waiver-request-status-icon
{
    height: 14px;
    width: 14px;
    margin: 0;
    padding: 0;
}

div.greenBar
{
    background-color: green;
}
div.amberBar
{
    background-color: orange;
}
div.redBar
{
    background-color: red;
}
div.trafficLight
{
    height: .75em;
    width: 3em;
    float: left;
    margin: .25em;
    padding: 0;
}
span.greenBar
{
    height: .75em;
    display: inline-block;
    background-color: green;
}
span.amberBar
{
    height: .75em;
    display: inline-block;
    background-color: orange;
}
span.redBar
{
    height: .75em;
    display: inline-block;
    background-color: red;
}

td.greenBar
{
    background-color: green !Important;
    width:20px;
    background-color: red;
}
td.amberBar
{
    background-color: orange !Important;
    width:20px;
}
td.redBar
{
    background-color: red !Important;
    width:20px;
}

.visible
{
    display: inherit;
}

.visible-inline
{
    display: inline;
}
.invisible
{
    display: none;
}

.badge-div
{
    display: inline-block;
    position: relative;
    line-height: normal;
    padding-right: 4px;
}
.badge
{
    background-color: #0000FF;
    color: white;
    font-weight: bold;
    font-size: 10px;
    width: 100%;
    height: 12px;
    text-align: center;
    display: inline-block;
    padding-top: 0px;
    padding-bottom: 2px;
    padding-left: 2px;
    padding-right: 2px;
    margin-left: 0;
    margin-right: 0;
    margin-top: 0;
    margin-bottom: 2px;
    border-radius: 45%;
}
.manageRUL-lnk a 
{
    vertical-align: top;
    height: 16px;
    width: 16px;
    background-image: url(./img/settings.png);
    background-repeat: no-repeat;
    background-position: right top;
    background-size: contain;
    display: inline-block;
    float: right;
}

.fakeLink
{
    color: blue;
}
.fakeLink:hover
{
    cursor: pointer;
}
#StatReviewerErrorMessage, #StatReviewerProcessingMessage
{
    display: none;
}
#StatReviewerErrorMessage tr td
{
    border: none;
}
}
