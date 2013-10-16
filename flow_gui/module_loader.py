import inspect
import re
from PyQt4 import QtGui, QtCore

from widgets import *

class GenericModuleWidget(QtGui.QWidget):

  state_changed = QtCore.pyqtSignal(int)
  
  def __init__(self, module):
    super(GenericModuleWidget, self).__init__()
    QtGui.QVBoxLayout(self)

    self.content = {}
    self.input_widgets_list = []

    self.element_selection_widget = QtGui.QComboBox()
    self.element_selection_widget.currentIndexChanged.connect(self.refresh)
    self.layout().addWidget(self.element_selection_widget)

    self.filter_module(module)
    self.module_name = module.__name__.rsplit('.',1)[-1]

    self.refresh()

  def filter_module(self, module):
    for factory_name, factory_object in inspect.getmembers(module):
      if factory_name[0]=='_' or re.search('Abstract', factory_name): # better way to do this with property?
        continue
      if inspect.isclass(factory_object) or inspect.isfunction(factory_object):
        self.append_factory(factory_name, factory_object)

  def append_factory(self, factory_name, factory_object):
    self.element_selection_widget.addItem(factory_name)
    input_widget = GenericIOWidget(factory_object)
    input_widget.output_generated.connect(self.set_contents)
    input_widget.input_changed.connect(self.check_state)

    self.input_widgets_list.append(input_widget)
    self.layout().addWidget(input_widget)

    self.refresh()

  def set_contents(self, current_input, output):
    self.last_input = current_input
    if isinstance(output, dict):
      self.content = output
    else:
      self.content = {self.module_name:output}
    self.check_state(current_input)

  def check_state(self, current_input):

    if not self.content:
      self.state_changed.emit(0)

    elif self.last_input == current_input:
      self.state_changed.emit(2)

    else:
      self.state_changed.emit(1)

  def refresh(self):
    current_index = self.element_selection_widget.currentIndex()

    for index, widget in enumerate(self.input_widgets_list):
      if index==current_index:
        self.visible_widget = widget
        self.visible_widget.show()

      else:
        widget.hide()