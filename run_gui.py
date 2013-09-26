import flow_gui
from OpenPNM import GUI

flow_gui.add_source(GUI.Network)
flow_gui.add_source(GUI.Algorithms)
flow_gui.run()