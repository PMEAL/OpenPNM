import OpenPNM
import flow_gui


flow_gui.add_source(OpenPNM.GUI.network)
flow_gui.add_source(OpenPNM.GUI.physics)
flow_gui.add_source(OpenPNM.GUI.algorithm)
flow_gui.add_source(OpenPNM.GUI.output)
flow_gui.run()