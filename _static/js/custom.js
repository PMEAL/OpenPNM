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

    // List of methods exclusive to Python's dict and list
    exclude = [
        'append', 'clear', 'copy', 'count', 'extend', 'fromkeys', 'get',
        'index', 'insert', 'items', 'keys', 'pop', 'popitem', 'remove',
        'reverse', 'setdefault', 'sort', 'update', 'values'
    ]
    // List of methods exclusive to numpy's ndaray
    exclude = exclude.concat(
        ['all', 'any', 'argmax', 'argmin', 'argpartition', 'argsort',
        'astype', 'base', 'byteswap', 'choose', 'clip', 'compress', 'conj',
        'conjugate', 'copy', 'ctypes', 'cumprod', 'cumsum', 'data',
        'diagonal', 'dot', 'dtype', 'dump', 'dumps', 'fill', 'flags', 'flat',
        'flatten', 'getfield', 'imag', 'item', 'itemset', 'itemsize', 'max',
        'mean', 'min', 'nbytes', 'ndim', 'newbyteorder', 'nonzero',
        'partition', 'prod', 'ptp', 'put', 'ravel', 'real', 'repeat',
        'reshape', 'resize', 'round', 'searchsorted', 'setfield', 'setflags',
        'shape', 'size', 'sort', 'squeeze', 'std', 'strides', 'sum',
        'swapaxes', 'take', 'tobytes', 'tofile', 'tolist', 'tostring',
        'trace', 'transpose', 'var', 'view', '__call__']
    )

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

    // Hide attributes exclusive to numpy's ndarray from tables
    exclude_attrs = [
        'base', 'ctypes', 'data', 'flags', 'flat', 'imag', 'itemsize',
        'nbytes', 'real', 'strides'
    ]
    for (let i = 0; i < exclude_attrs.length; i++) {
        // Search through toctree for "exclude_attrs[i]"
        tmp = $(`dt:contains(${exclude_attrs[i]})`)
        // Hide that element
        tmp.hide()
        // Hide attr's description (2nd column)
        tmp.next("dd").hide()
    }

    // Change h1 font to monospace if inside Module Reference
    tmp = $("a:contains('Module Reference')")[0]
    if (tmp.parentElement.classList.contains("active") == true){
        $("h1").css("font-family", "Ubuntu Mono")
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
