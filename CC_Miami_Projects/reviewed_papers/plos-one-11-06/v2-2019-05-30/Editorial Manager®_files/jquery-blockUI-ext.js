//----------------------------------------------------------------------------------------
//  Copyright © 2013-Present Aries Systems Corporation. All Rights Reserved.
//  Copying, reverse engineering, adaptation or any other derivative use
//  prohibited. This material is proprietary and confidential information
//  of Aries Systems Corporation.
//  
//  Date Created: 20130613
//  Version Introduced: 10.0
//
//  Description:  Defaults and variables for block UI plug-in  
//
//----------------------------------------------------------------------------------------
//
//-------------- Spec 13.0-37 --------------
// 20151108  CSR
// Replaced contents of this file with code found in jquery-blockUI-ext2.js. In addition, defunct 
// file jquery-blockUI-ext2.js to remove equivalent functionality for maintainability purposes.
// 
// -------------- Bugfix 26493 --------------
// 20131224  GBS
// Added blockUISpinner function to output overlay with indictor image.
// Added spinnerCSS for formatting.
// Created $.blockUISpinner.spinnerimg to hold preloaded image based on location
//

(function() {
    // various options for $.blockUI plugin
    var options = {
        // just plain gray overlay
        plain: {
            message: '',
            overlayCSS: {
                backgroundColor: '#000',
                opacity: 0.2,
                cursor: 'wait'
            }
        },
        // spinner without transparent overlay
        spinnerTransparentOverlay: {
            message: '<div class="spinner"></div>',
            css: {
                border: 'none',
                top: '25%',
                backgroundColor: 'transparent',
                cursor: 'default'
            },
            overlayCSS: {
                backgroundColor: '#000',
                opacity: 0.0,
                cursor: 'default'
            },
            messageDelay: 1000,
            showOverlay: true
        },
        // spinner with overlay and message text
        spinnerOverlayWithMessage: {
            message: '<div class="spinner-wrapper"><div class="spinner"></div><div class="spinner-message">{0}</div></div>',
            css: {
                border: 'none',
                top: '25%',
                backgroundColor: 'transparent',
                cursor: 'default'
            },
            centerX: false,
            cetnerY: false,
            overlayCSS: {
                backgroundColor: '#000',
                opacity: 0.2,
                cursor: 'default'
            },
            showOverlay: true
        }
    };

    // wrapper functions
    // compatible options can be combined as: $.blockUI2('option1 option2')
    $.blockUI2 = function(opts, text) {
        if (typeof opts === 'string') {
            var optStrs = opts.split(/\s+/),
                m = {};

            $.each(optStrs, function (i, s) {
                // insert message text for spinner with message
                if (s == 'spinnerOverlayWithMessage') {
                    options[s].message = options[s].message.replace('{0}', text);
                }
                m = $.extend(true, m, options[s] || {});
            });

            // not valid option(s), fall back to 'plain'
            if ($.isEmptyObject(m)) {
                m = options['plain'];
            }
            opts = m;
        }
        if (!opts) {
            opts = options['plain'];
        }
        $.blockUI(opts);
    };

    $.unblockUI2 = function(opts) {
        $.unblockUI(opts);
    };
})();
