from __future__ import absolute_import, print_function

import inspect
import logging

from PyQt4 import QtGui, QtCore

import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

def required_arguments(callable_object):
  arg_names, args_packer, kwargs_packer, default_values = inspect.getargspec(callable_object)
  if not default_values: default_values = []
  return arg_names[:-len(default_values)]


class MatplotlibWidget(FigureCanvas):

  def __init__(self, vis_method):
    super(MatplotlibWidget, self).__init__(Figure())
    self.vis_method = vis_method

  def plot(self, available_data):
    self.figure.clf()
    input_data = dict( (key, available_data[key]) for key in required_arguments(self.vis_method) )
    logging.debug("Obtained these input arguments: {0}".format(input_data) )    
    self.vis_method(fig=self.figure, **input_data)


class PreviewWidget(QtGui.QTabWidget):
  
  def __init__(self, preview_source=None):
    super(PreviewWidget, self).__init__()

    self.tab_widgets = []
    for vis_name, vis_method in inspect.getmembers(preview_source, inspect.isfunction):
      self.tab_widgets.append( (vis_name, MatplotlibWidget(vis_method)) )

    self.debug_tab = QtGui.QTextEdit()
    self.debug_tab.setAcceptRichText(False)
    self.debug_tab.setFocusPolicy(QtCore.Qt.NoFocus)
    self.addTab(self.debug_tab, "Debug")

  def update(self, property_dict):
    ''' go through all of the supplied preview sources
    
    if the arguments are available, supply all relevant visuals
    '''
    for (name, w) in self.tab_widgets:
      idx = self.indexOf(w)
      try:
        w.plot(property_dict)
        self.insertTab(0, w, name)

      except Exception as e:
        if idx >= 0:
          self.removeTab(idx)
        logging.debug("Attempted to plot {name} and failed because {e}".format(**locals()))

    self.debug_tab.setPlainText(
      '\n\n'.join("{0}:\n{1}".format(key, value) for key, value in sorted(property_dict.items())) )
    


if __name__ == '__main__':
  import sys; sys.path.append(__file__.rsplit('/',2)[0])
  import run_gui