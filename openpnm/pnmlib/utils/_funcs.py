import time
import sys
import numpy as np


__all__ = [
    'tic',
    'toc',
]


def _format_time(timespan, precision=3):
    """Formats the timespan in a human readable form"""

    if timespan >= 60.0:
        # we have more than a minute, format that in a human readable form
        # Idea from http://snipplr.com/view/5713/
        parts = [("d", 60*60*24), ("h", 60*60), ("min", 60), ("s", 1)]
        time = []
        leftover = timespan
        for suffix, length in parts:
            value = int(leftover / length)
            if value > 0:
                leftover = leftover % length
                time.append(u'%s%s' % (str(value), suffix))
            if leftover < 1:
                break
        return " ".join(time)

    # Unfortunately the unicode 'micro' symbol can cause problems in
    # certain terminals.
    # See bug: https://bugs.launchpad.net/ipython/+bug/348466
    # Try to prevent crashes by being more secure than it needs to
    # E.g. eclipse is able to print a µ, but has no sys.stdout.encoding set.
    units = [u"s", u"ms", u'us', "ns"]  # the save value
    if hasattr(sys.stdout, 'encoding') and sys.stdout.encoding:
        try:
            u'\xb5'.encode(sys.stdout.encoding)
            units = [u"s", u"ms", u'\xb5s', "ns"]
        except UnicodeEncodeError:
            pass
    scaling = [1, 1e3, 1e6, 1e9]

    if timespan > 0.0:
        order = min(-int(np.floor(np.log10(timespan)) // 3), 3)
    else:
        order = 3
    return u"%.*g %s" % (precision, timespan * scaling[order], units[order])


def tic():
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    See Also
    --------
    toc

    """
    global _startTime_for_tictoc
    _startTime_for_tictoc = time.time()


def toc(quiet=False):
    r"""
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    Parameters
    ----------
    quiet : bool, default is False
        If False then a message is output to the console. If
        True the message is not displayed and the elapsed time is returned.

    See Also
    --------
    tic

    """
    if "_startTime_for_tictoc" not in globals():
        raise Exception("Start time not set, call tic first")
    t = time.time() - _startTime_for_tictoc
    if quiet is False:
        print(f"Elapsed time: {_format_time(t)}")
    return t
