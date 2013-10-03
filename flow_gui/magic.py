import inspect
import re
from decimal import Decimal
from PyQt4 import QtGui, QtCore

class GenericIOWidget(QtGui.QWidget):
  output_generated = QtCore.pyqtSignal(dict, object)
  input_changed = QtCore.pyqtSignal(dict)

  def __init__(self, factory):
    self.factory = factory

    super(GenericIOWidget, self).__init__()
    QtGui.QVBoxLayout(self)

    self.input_form = QtGui.QFormLayout()
    self.layout().addLayout(self.input_form)

    self.create_input_widgets()

    output_button = QtGui.QPushButton("Process")
    output_button.clicked.connect(self.current_output)
    self.layout().addWidget(output_button)

  def create_input_widgets(self):
    # create persistent container for reference
    self.widget_list = []

    # determine where to extract the arguments from
    if inspect.isclass(self.factory):
      self.arg_names_list, arg_values_list = inspect.getargspec(self.factory.__init__)[::3]
      self.arg_names_list.remove('self')

    elif inspect.isfunction(self.factory):
      self.arg_names_list, arg_values_list = inspect.getargspec(self.factory)[::3]
    
    # determine the number of state inputs. account for possible empty list
    if arg_values_list is None: arg_values_list = []
    n_state_inputs = len(self.arg_names_list) - len(arg_values_list)
    
    # state-input variables
    for arg_name in self.arg_names_list[:n_state_inputs]:
      widget = QtGui.QWidget()
      widget.get = lambda arg_name=arg_name: self.parent().item.branch_properties().get(arg_name)
      self.widget_list.append(widget)
    
    # editable variables
    for arg_name, arg_value in zip(self.arg_names_list[n_state_inputs:], arg_values_list):
      widget = self.create_widget(arg_name, arg_value)
      self.widget_list.append(widget)
      label = arg_name
      self.input_form.addRow(label, widget)

  def create_widget(self, arg_name, arg_value=None):

    if isinstance(arg_value, bool):
      widget = QtGui.QCheckBox()
      widget.setCheckState(arg_value)
      widget.get = lambda: bool(widget.checkState())

      widget.stateChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, int):
      widget = QtGui.QSpinBox()
      widget.setRange(0,1000)
      widget.setValue(arg_value)
      widget.get = widget.value

      widget.valueChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, (float)):
      widget = QtGui.QDoubleSpinBox()
      exp = Decimal(str(arg_value)).as_tuple().exponent
      widget.setDecimals(-exp)
      widget.setSingleStep(10**exp)
      widget.setValue(arg_value)
      widget.get = widget.value

      widget.valueChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif arg_name=='path':
      widget = QtGui.QWidget()
      hbox = QtGui.QHBoxLayout(widget)
      line_edit = QtGui.QLineEdit()
      line_edit.setReadOnly(True)
      line_edit.setText(arg_value)
      hbox.addWidget(line_edit)
      find_button = QtGui.QPushButton()
      find_button.clicked.connect(lambda: line_edit.setText(QtGui.QFileDialog.getOpenFileName()[0]))
      hbox.addWidget(find_button)
      widget.get = lambda: line_edit.text()

      line_edit.textChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, (str)):
      widget = QtGui.QLineEdit()
      widget.setText(arg_value)
      widget.get = lambda: widget.text()

      widget.textEdited.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, GenericModuleWidget):
      widget = QtGui.QComboBox()
      widget.setModel(arg_value.item_model)
      widget.get = lambda: arg_value.item_model.item(widget.currentIndex()).data().toPyObject()

      widget.currentIndexChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif hasattr(arg_value, "__iter__"):
      widget = QtGui.QWidget()
      widget.setLayout(QtGui.QHBoxLayout())
      widget.sub_widgets = [self.create_widget(None, sub_value) for sub_value in arg_value]
      widget.get = lambda: [sub_widget.get() for sub_widget in widget.sub_widgets]

      for sub_widget in widget.sub_widgets:
        widget.layout().addWidget(sub_widget)

    else:
      widget = QtGui.QLabel("%s." % type(arg_value))
      widget.get = lambda: None

    return widget

  def current_input(self):
    arg_values = [widget.get() for widget in self.widget_list]
    return dict(zip(self.arg_names_list, arg_values))

  def current_output(self):
    self.output_generated.emit(
      self.current_input(),
      self.factory( **self.current_input() )
      )

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

# TEST STUFF GOES HERE
if __name__ == '__main__':
  import sys; sys.path.append(__file__.rsplit('/',2)[0]+'/modules')
  import simulation

  app = QtGui.QApplication([])

  pg = GenericModuleWidget(simulation)
  pg.show()

  app.exec_()