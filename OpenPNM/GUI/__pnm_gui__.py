# -*- coding: utf-8 -*-

r"""
*************************************************************************
:mod:`OpenPNM.GUI.__pnm_gui__` -- Graphical User Interfaces (GUI)
*************************************************************************

.. module:: OpenPNM.GUI.__pnm_gui__

Contents
========
The OpenPNM package imports all the functions from the top level modules. 
 
Import
======
>>> import OpenPNM as PNM
>>> PNM.GUI.run_pnm_gui()


Submodules
==========


Classes
=======

"""


from PyQt4 import QtGui, QtCore

TITLE = "Prototype"

class BrowseFileWidget(QtGui.QWidget):
  def __init__(self, label=None):
    super(BrowseFileWidget,self).__init__()
    hbox = QtGui.QHBoxLayout(self)

    if label:
      hbox.addWidget(QtGui.QLabel(label))

    self.path_widget = QtGui.QLineEdit()
    hbox.addWidget(self.path_widget)
    self.path_widget.setReadOnly(True)
    browse_button = QtGui.QPushButton("Browse")
    browse_button.clicked.connect(self.browse)
    hbox.addWidget(browse_button)

  def browse(self):
    path = QtGui.QFileDialog.getOpenFileName()

    if path:
      self.path_widget.setText(path)

class PropertySpinBox(QtGui.QWidget):
  def __init__(self, label, default, limits):
    super(PropertySpinBox,self).__init__()
    hbox = QtGui.QHBoxLayout(self)

    if isinstance(default,int):
      self.spinbox = QtGui.QSpinBox()
    else:
      self.spinbox = QtGui.QDoubleSpinBox()

    self.spinbox.setValue(default)
    hbox.addWidget(self.spinbox)
    hbox.addWidget(QtGui.QLabel(label))

  def value(self):
    return self.spinbox.value()

class CubicSettingsPanel(QtGui.QWidget):
  def __init__(self):
    super(CubicSettingsPanel,self).__init__()
    self.settings = {}
    grid = QtGui.QGridLayout(self)

    for (label, default,  limits,   row, col, key) in [
        ('X',    10,      [1,100],  0,   0  , 'x'),
        ('Y',    10,      [1,100],  0,   1  , 'y'),
        ('Z',    10,      [1,100],  0,   2  , 'z'),
        ('alpha',    3.,      [0,100],  1,   0  , 'wshape'),
        ('beta',    1.,      [0,100],  1,   1  , 'wscale'),
        ('c/o',  4.,      [0,100],  1,   2  , 'cut-off'),]:

        self.settings[key] = PropertySpinBox(label, default, limits)
        grid.addWidget(self.settings[key], row, col)

  def get(self):
    return self.settings






class GenericTab(QtGui.QWidget):
  def __init__(self):
    super(GenericTab,self).__init__()

class NetworkTab(GenericTab):

  model_generated = QtCore.pyqtSignal(str)

  def __init__(self, parent):
    super(NetworkTab,self).__init__()
    parent.addTab(self, "Network")

    grid = QtGui.QGridLayout(self)

    grid.addWidget(QtGui.QLabel("Name"),0,0)
    self.name_widget = QtGui.QLineEdit("New Network")
    grid.addWidget(self.name_widget,0,1,1,-1)

    self.import_radio_button = QtGui.QRadioButton("Import")
    grid.addWidget(self.import_radio_button,1,0)
    self.import_path_widget = BrowseFileWidget()
    grid.addWidget(self.import_path_widget,1,1)

    cubic_radio_button = QtGui.QRadioButton("3D Cubic")
    cubic_radio_button.setChecked(True)
    grid.addWidget(cubic_radio_button,2,0)
    cubic_settings_panel = CubicSettingsPanel()
    grid.addWidget(cubic_settings_panel,2,1,1,-1)

    preview_plot = QtGui.QLabel("Plot")
    preview_plot.setAlignment(QtCore.Qt.AlignCenter)
    preview_plot.setMinimumSize(240,180)
    grid.addWidget(preview_plot,5,0,1,-1)

    generate_button = QtGui.QPushButton("Generate Model")
    generate_button.clicked.connect(self.model_generate)
    grid.addWidget(generate_button,6,0,1,-1)

  def model_generate(self):
    self.model_generated.emit(self.name_widget.text())

class SimulationTab(GenericTab):
  def __init__(self, parent):
    super(SimulationTab,self).__init__()
    parent.addTab(self, "Simulation")

    vbox = QtGui.QVBoxLayout(self)

    inlet_box = QtGui.QGroupBox("Inlet")
    inlet_box.setLayout(QtGui.QVBoxLayout())
    inlet_box.layout().addWidget(QtGui.QRadioButton("full-face [flux]"))
    inlet_box.layout().addWidget(QtGui.QRadioButton("full-face [pressure]"))
    inlet_box.layout().addWidget(QtGui.QRadioButton("point"))
    vbox.addWidget(inlet_box)
    
    alg_box = QtGui.QGroupBox("Algorithm")
    alg_box.setLayout(QtGui.QHBoxLayout())
    alg_box.layout().addWidget(QtGui.QRadioButton("invasion\npercolation"))
    alg_box.layout().addWidget(QtGui.QRadioButton("access-limited\nordinary percolation"))
    vbox.addWidget(alg_box)
    
    self.capillary_file_widget = BrowseFileWidget("Capillary Pressure File")
    vbox.addWidget(self.capillary_file_widget)

class TransportTab(GenericTab):
  def __init__(self, parent):
    super(TransportTab,self).__init__()
    parent.addTab(self, "Transport")

    form = QtGui.QFormLayout(self)
    form.setVerticalSpacing(20)

    for mechanism in ['Diffusion', 'Permeability', 'Through-plane', 'In-plane']:
      form.addRow(QtGui.QCheckBox(mechanism))

class VisualizationTab(GenericTab):
  def __init__(self, parent):
    super(VisualizationTab,self).__init__()
    parent.addTab(self, "Visualization")

    form = QtGui.QFormLayout(self)
    form.setVerticalSpacing(20)

    mechanism_box = QtGui.QComboBox()
    mechanism_box.addItems(['Concentration', 'Pressure', 'Filling pressure', 'Filling order'])
    form.addRow("Color", mechanism_box)

    preview_plot = QtGui.QLabel("3D Plot")
    preview_plot.setAlignment(QtCore.Qt.AlignCenter)
    preview_plot.setMinimumSize(240,240)
    form.addRow(preview_plot)

class MainWindow(QtGui.QMainWindow):

  def __init__(self):
    super(MainWindow,self).__init__()

    self.setWindowTitle(TITLE)
    
    self.tab_widget = QtGui.QTabWidget(self)
    self.setCentralWidget(self.tab_widget)

    self.net_tab = NetworkTab(self.tab_widget)
    self.sim_tab = SimulationTab(self.tab_widget)
    self.tp_tab = TransportTab(self.tab_widget)
    self.vis_tab = VisualizationTab(self.tab_widget)

    self.model_reference = QtGui.QLineEdit()
    self.model_reference.setReadOnly(True)
    self.toolbar = QtGui.QToolBar()
    self.toolbar.addWidget(self.model_reference)
    self.addToolBar(self.toolbar)

    self.net_tab.model_generated.connect(self.model_reference.setText)
    
def run_pnm_gui():
  import sys

  app = QtGui.QApplication(sys.argv)

  main = MainWindow()
  main.show()

  sys.exit(app.exec_())

if __name__ == '__main__':
  run_pnm_gui()

