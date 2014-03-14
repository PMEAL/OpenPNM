 .. _sandbox: 'Table of contents for Sandbox'


###############################################################
Sandbox - Use pound symbol (above and below) to indicate a part
###############################################################

++++++++++++++++++++++++++++++++++++++++++++++++++++++
Chapters are identified by asterisks (above and below)
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Sections are identified by equal signs placed underneath
========================================================

Subsections are identified by dashes placed underneath
------------------------------------------------------

To give an example expression, we follow text by two colons.  For example::

  .. math::

    s_{\textrm{wp}}=s^{*}_{\textrm{wp}} \left (\frac{P^{*}_{\textrm{c}}}{P_{\textrm{c}}} \right )^{\eta}, \ P_{\textrm{c}}>P^{*}_{\textrm{c}}

A pdf for writing Latex math expressions:
ftp://ftp.ams.org/pub/tex/doc/amsmath/amsldoc.pdf

A very handy link for testing Latex math expressions:
http://www.codecogs.com/latex/eqneditor.php

+++++++++++++++++
Late pore filling
+++++++++++++++++

Due to the presence of acute internal features, the pressure required to completely fill 
a pore with the non-wetting phase may be higher than the pressure required for entry. The 
following expression is used to model the late pore filling phenomenon:

.. math::

   s_{\text{wp}}=s^{*}_{\text{wp}} \left (\frac{P^{*}_{\text{c}}}{P_{\text{c}}} \right )^{\eta}, \ P_{\text{c}}>P^{*}_{\text{c}}

   
.. math::

   s_{\textrm{wp}}
   
where :math:`\eta` is the filling exponent, :math:`s_{\text{wp}}` is the wetting phase 
saturation of a given pore at capillary pressure :math:`P_{\text{c}}`, and :math:`s^{*}_{\text{wp}}` 
is the wetting phase saturation of the same pore at the capillary pressure, 
:math:`P^{*}_{\text{c}}`, corresponding to first entry (breakthrough) of the non-wetting phase. 
The parameters :math:`{\eta}` and :math:`s^{*}_{\text{wp}}` are adjustable.