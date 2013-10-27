import OpenPNM

import flow_gui

flow_gui.add_source(OpenPNM.GUI.network)
flow_gui.add_source(OpenPNM.GUI.fluid)
# flow_gui.add_source(OpenPNM.Physics.CapillaryPressure, name="physics")
flow_gui.add_source(OpenPNM.GUI.algorithm)
# flow_gui.add_source(OpenPNM.GUI.output)
flow_gui.set_previews(OpenPNM.Visualization.Plots)
flow_gui.run()