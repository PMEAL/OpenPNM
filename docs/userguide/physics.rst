. _physics:

###############################################################################
Pore Scale Physics
###############################################################################

In OpenPNM each physical process is given its own module, which are defined by the files stored in /OpenPNM/Physics.  Within each module any number of model equations can be defined by simply adding the appropriate method/equations.  This setup was designed to facilitate easy editing of the physics models which can be highly varied, and are likely to be the most customized aspect of a given researchers work.  

For a complete listing of the existing modules and their included methods see the :ref:`function_ref`.

===============================================================================
Adding Custom Models to Existing Physics Classes
===============================================================================
Defining custom models for different physical processes is designed to be as easy as possible in OpenPNM. There are, however, a few technical and convention considerations to bear in mind. The fist two are relatively obvious, but the others are very important.

1. All methods should take the Fluid Object to which the physics are being applied as well as the Network Object from which physical structure is taken.  

2. Methods should 'know' what value(s) they are calculating and they should add or update the necessary pore and/or throat conditions directly on the Fluid Object.  

3. In general, when methods require physical conditions in the Fluid (such as temperature to calculate viscosity) these conditions must already be assigned to the fluid.  That is, before calculating a temperature dependent viscosity it stands to reason that the temperature should be set.  (Note that when a Fluid is created it is given a standard conditions by default)

4. Methods must be able to handle both scalar and vector arguments.  Numpy's broadcasting abilities make this almost transparent, with errors occurring only when vectors of different length are broadcast together.  

.. warning:: 
    Although scalars and vectors are broadcast together by Numpy without issue, problems can arise if a user attempts to use indexing (like in a for loop).  This will work fine on vectors, but will fail for scalars.  Neither Python or Numpy will index into a 0 dimensional variable (i.e. a scalar).  

Consider the example of calculating capillary entry pressures with the Washburn equation.  The call to this function looks like:

.. code-block:: python
    
    OpenPNM.Physics.CapillaryPressure.Washburn(pn, water)
    
And the method itself looks like: 

.. code-block:: python
    
	sigma = water.throat_conditions['surface_tension']
	theta = water.throat_conditions['contact_angle']
    Pc = -4*sigma*sp.cos(sp.radians(theta))/net.throat_properties['diameter']
    
If the network is isothermal then sigma will be constant throughout, and only a scalar value need be sent.  If the network is nonisothermal, then sigma will vary with the local temperature.  Assuming that a value of sigma as a function of temperature has been calculated for each throat and stored in pn.throat_properties['sigma'], then the above call would look like:

.. code-block:: python
    
    OpenPNM.Physics.CapillaryPressure.Washburn(pn, water)
    
where theta is still a scalar.  The equation for finding Pc is completely unchanged.  The equation was already vectorized due to the throat diameter in the demoninator.  The presence of a vector in the numerator will be handled transparently and efficiently by Numpy.  Importantly, OpenPNM allows for storage of scalar values in the pore and throat dictionaries on the assumption that this means these values apply to all pores and throats in the network.  This means that in an isothermal network sigma can be stored as a scalar throat property and the Washburn equation called using the second option above.  As can be seen, this means the code is completely indifferent to scalar or vector arguments.  

And finally, it is highly encouraged that quality documentation is added (see next section).

-------------------------------------------------------------------------------
Documenting Your Method
-------------------------------------------------------------------------------
OpenPNM uses the Numpydoc standard of documentation.  This standard defines what information should be included in the documentation of each method, class and module.  The advantage of Numpydoc is the vast majority of numerical and array function used by OpenPNM are based on Numpy.  By keeping the documentation style consistent it provides a more seamless work flow.  For instance when working in certain IDEs (such as Spyder), the 'Object Inspector' shows the docstring of method currently being entered on the console.  By using the Numpydoc style, the information in the Object Inspector is consistent between Numpy (Scipy) and OpenPNM.  

The second aspect to consider when documenting a physics module is that the Numpy 'Notes' section is used heaving in OpenPNM.  In fact, the aim is to have no model-specific documentation written in these pages; instead it will all be automatically generated by pulling the docstring from the code.  All the information above was generated this way.  

===============================================================================
Adding New Physics Classes to the Module
===============================================================================
If the supplied physics modules are not sufficient for some esoteric type of simulation, it is not difficult add additional modules.  This is accomplished in two steps.  First, a file with the desired physics name is added to the /OpenPNM/Physics folder.  Naturally, this file will be populated with the desired methods.  There is no special requirements on the contents of this file, such as class definitions or __init__ methods.  There are, however, a few steps required to register this new module with OpenPNM.  In the file /OpenPNM/Physics/__init__.py, the module must be imported.  If the module name is Acoustics for instance, then the following must be added to the __init__.py file:

.. code-block:: python
    
    import acoustics
    
This addition means that the Acoustics module and the methods therein will be under the OpenPNM.Physics namespace.  

To ensure this module is properly included in the documentation is a bit more convoluted.  



