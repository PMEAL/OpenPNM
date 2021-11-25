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

    // var selected = $("a.current")
    // selected.html("&#8594; " + selected.text())

    // $(".toctree-l3").hide()

    exclude = [
        'append', 'clear', 'copy', 'count', 'extend', 'fromkeys', 'get',
        'index', 'insert', 'items', 'keys', 'pop', 'popitem', 'remove',
        'reverse', 'setdefault', 'sort', 'update', 'values'
    ]

    // Hide methods exclusive to dict and list from tables
    for (let i = 0; i < exclude.length; i++) {
        // Search through first column of the table for "exclude[i]("
        tmp = $("tr").find(`td:first:contains(${exclude[i]}()`)
        // Find the row(s) containing the returned query
        row = tmp.closest('tr')
        // Hide that row
        row.hide()
        // Comment the line above and uncomment the next for DEBUGGING
        // row.css("background-color", "red")
    }

    // Hide methods exclusive to dict and list from toctree
    for (let i = 0; i < exclude.length; i++) {
        // Search through toctree for "exclude[i]"
        tmp = $(`.toctree-l3:contains(${exclude[i]})`)
        // Hide that row
        tmp.hide()
        // Comment the line above and uncomment the next for DEBUGGING
        // tmp.css("background-color", "red")
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
