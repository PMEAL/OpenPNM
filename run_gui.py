import flow_gui
from OpenPNM import GUI

flow_gui.add_source(GUI.network)
flow_gui.add_source(GUI.algorithm)
flow_gui.run()