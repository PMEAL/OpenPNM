// We don't want the sidebar to include list of examples, that's just
// much. This simple script tries to select those items and hide them.
// "aside" selects sidebar elements, and "href" narrows it down to the
// list of examples. This is a workaround, not a permanent fix.

// The next lines work with furo theme
// var examples = $('aside a[href*="examples.html"]')
// var examples_clicked = $( ":contains('Examples')" ).filter($( ".current.reference.internal" ))
// examples.nextAll().hide()
// examples_clicked.nextAll().hide()

// The next lines work with pydata theme

if (location.protocol.startsWith("http") & location.protocol !== 'https:') {
    location.replace(`https:${location.href.substring(location.protocol.length)}`);
}

window.onload = function () {
    var examples_clicked = $( ".active a:contains(Examples)" )
    if (examples_clicked.length == 1) {
        $(" nav.bd-links ").children().hide()
    }
};

// window.onload = function() {
//     if (window.jQuery) {
//         // jQuery is loaded
//         alert("Yeah!");
//     } else {
//         // jQuery is not loaded
//         alert("Doesn't Work");
//     }
// }
