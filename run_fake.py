import fake

import flow_gui
flow_gui.add_source(fake.network)
flow_gui.add_source(fake.thermophysics, name='fluid')
flow_gui.add_source(fake.physics)
# flow_gui.add_source(OpenPNM.Physics.CapillaryPressure, name="physics")
# flow_gui.add_source(OpenPNM.GUI.algorithm)
# flow_gui.add_source(OpenPNM.GUI.output)
flow_gui.run()