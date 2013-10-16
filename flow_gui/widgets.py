from __future__ import print_function

import inspect
from decimal import Decimal

from PyQt4 import QtGui, QtCore


class GenericIOWidget(QtGui.QWidget):
  output_generated = QtCore.pyqtSignal(dict, object)
  input_changed = QtCore.pyqtSignal(dict)

  def __init__(self, factory, button=True):
    self.factory = factory

    super(GenericIOWidget, self).__init__()
    QtGui.QVBoxLayout(self)

    self.input_form = QtGui.QFormLayout()
    self.layout().addLayout(self.input_form)

    self.create_input_widgets()

    if button:
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

    elif isinstance(self.factory, dict):
      self.arg_names_list, arg_values_list = zip(*sorted(self.factory.items()))
    
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

    elif arg_name.split('_')[0]=='stats':
      widget = SciPyStatsWidget()
      widget.changed.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, (str)):
      widget = QtGui.QLineEdit()
      widget.setText(arg_value)
      widget.get = lambda: widget.text()

      widget.textEdited.connect(lambda: self.input_changed.emit(self.current_input()) )

    # elif isinstance(arg_value, GenericModuleWidget):
    #   widget = QtGui.QComboBox()
    #   widget.setModel(arg_value.item_model)
    #   widget.get = lambda: arg_value.item_model.item(widget.currentIndex()).data().toPyObject()

    #   widget.currentIndexChanged.connect(lambda: self.input_changed.emit(self.current_input()) )

    elif isinstance(arg_value, (list, tuple)):
      ''' if the input is iterable, attempt to capture its structure by creating a widget nest
      as we travel through the tree
      '''
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


class SciPyStatsWidget(QtGui.QWidget):
  changed = QtCore.pyqtSignal()
  
  def __init__(self):
    super(SciPyStatsWidget, self).__init__()
    QtGui.QVBoxLayout(self)

    self.distribution_name = QtGui.QComboBox()
    self.distribution_name.addItems('weibull_min gamma normal lognormal'.split())
    self.layout().addWidget(self.distribution_name)

    self.parameters = GenericIOWidget({'shape':1.5, 'loc':6E-6, 'scale':2E-5}, button=False)
    self.layout().addWidget(self.parameters)

  def get(self):
    return dict([('name', str(self.distribution_name.currentText()))] + self.parameters.current_input().items())


if __name__ == '__main__':
  app = QtGui.QApplication([])

  g = GenericIOWidget(lambda stats_some='something': stats_some)
  g.output_generated.connect(print)
  g.show()

  app.exec_()