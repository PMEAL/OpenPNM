import sip
from PyQt4 import QtGui, QtCore
from .module_loader import GenericModuleWidget
from . import icons_rc

VERSION = 0.1
NAME = "OpenPNM"

class PipelineItem(QtGui.QStandardItem):

  def __init__(self, text, module_widget, main_window):
    super(PipelineItem, self).__init__(text)
    self.setTristate(True)
    self.setCheckState(0)
    self.module_widget = module_widget
    self.module_widget.item = self
    
    self.module_widget.state_changed.connect(self.update)

    self.preview = QtGui.QTextEdit()
    self.preview.setAcceptRichText(False)
    self.preview.setFocusPolicy(QtCore.Qt.NoFocus)

    self.main_window = main_window

  def update(self, check_state):
    self.setCheckState(check_state)
    self.preview.setPlainText(
      '\n\n'.join("{0}:\n{1}".format(key, value) for key, value in sorted(self.branch_properties().items())) )
    self.main_window.update()

  @property
  def data_dict(self):
    return self.module_widget.content

  def branch_properties(self):
    item = self
    properties = {}

    while item:
      for key, value in item.data_dict.items():
        properties.setdefault(key, value)
      item = item.parent()

    # for key, value in DEFAULTS.items():
      # properties.setdefault(key, value)

    return properties

class MainWindow(QtGui.QMainWindow):

  def __init__(self):
    super(MainWindow, self).__init__()
    self.setWindowTitle("{name} v{version}".format(name=NAME,version=VERSION))
    self.setCorner( QtCore.Qt.BottomRightCorner, QtCore.Qt.RightDockWidgetArea )
    self.persistent_container = [] # avoid quirk with model items being auto-deleted silently by Qt
    self.create_model()
    self.create_widgets()
    self.create_actions()
    self.create_menubar()

  def create_model(self):
    self.model = QtGui.QStandardItemModel()

  def create_actions(self):
    self.remove_action        = QtGui.QAction("Remove", self, triggered=self.remove_branch, shortcut=QtGui.QKeySequence.Delete, icon=QtGui.QIcon(":qrc/icons/delete.png"))
    self.run_all_action       = QtGui.QAction("Run All", self, triggered=self.run_all, shortcut="Ctrl+R", icon=QtGui.QIcon(":qrc/icons/resultset_next.png"))
    self.save_project         = QtGui.QAction("Save Project", self, triggered=self.save_project, shortcut="Ctrl+S", icon=QtGui.QIcon(":qrc/icons/disk.png"))

    self.view_settings_action = self.settings.toggleViewAction()
    self.view_settings_action.setIcon(QtGui.QIcon(":qrc/icons/cog.png"))
    self.view_visualization_action = self.preview.toggleViewAction()
    self.view_visualization_action.setIcon(QtGui.QIcon(":qrc/icons/monitor.png"))

  def create_widgets(self):
    self.toolbar = QtGui.QToolBar(self)
    self.toolbar.setFloatable(False)
    self.toolbar.setMaximumHeight(32)
    self.addToolBar(self.toolbar)

    self.view = QtGui.QTreeView()
    self.view.setModel(self.model)
    self.view.setHeaderHidden(True)
    self.view.selectionModel().selectionChanged.connect(self.update)
    self.view.setMinimumWidth(400)
    self.view.setStyleSheet('''
      QTreeView::indicator                { width: 16px; height: 16px }
      QTreeView::indicator:checked        { background-image:url(":qrc/icons/tick.png") }
      QTreeView::indicator:unchecked      { background-image:url(":qrc/icons/cross.png") }
      QTreeView::indicator:indeterminate  { background-image:url(":qrc/icons/error.png") }
      '''
      )
    self.setCentralWidget(self.view)

    self.settings = QtGui.QDockWidget("Settings")
    self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.settings)

    self.preview = QtGui.QDockWidget("Visualization")
    self.preview.setMinimumSize(500,500)
    self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.preview)

  def populate_toolbar(self):
    ''' There is no way to add items to the beginning of the toolbar.
        This helper function allows us to append the standard actions
        as the window is shown, by calling it right before the `show()`
        method is called.
    '''
    self.toolbar.addAction(self.remove_action)
    self.toolbar.addAction(self.run_all_action)
    self.toolbar.addSeparator()
    self.toolbar.addAction(self.view_settings_action)
    self.toolbar.addAction(self.view_visualization_action)
    self.toolbar.addSeparator()
    self.toolbar.addAction(self.save_project)

  def show(self, *args, **kwargs):
    self.populate_toolbar()
    self.statusBar().showMessage("Ready.")
    super(MainWindow, self).show(*args, **kwargs)

  def create_menubar(self):
    self.menubar = self.menuBar()

    self.filemenu = self.menubar.addMenu("&File")
    self.editmenu = self.menubar.addMenu("&Edit")
    self.networkmenu = self.menubar.addMenu("&Network")
    self.simmenu = self.menubar.addMenu("&Simulation")
    self.helpmenu = self.menubar.addMenu("&Help")

  def add_source(self, module, name=None, icon=None, shortcut=None):
    if name:
      module_name = name
    else:
      module_name = module.__name__.rsplit('.')[-1]

    def anonymous_build_method():
      parent_index, parent_item = self.selection()
      if not parent_item:
        parent_item = self.model.invisibleRootItem()

      new_item = PipelineItem("New {type}".format(type=module_name),
                              GenericModuleWidget(module),
                              main_window=self)
      self.persistent_container.append(new_item)
      parent_item.appendRow(new_item)
      if parent_index:
        self.view.expand(parent_index)
      self.view.selectionModel().select(new_item.index(), QtGui.QItemSelectionModel.ClearAndSelect)

    new_action = QtGui.QAction("Add {name}".format(name=module_name), self, triggered=anonymous_build_method)
    if icon:      new_action.setIcon(QtGui.QIcon(icon))
    if shortcut:  new_action.setShortcut(shortcut)
    self.toolbar.addAction(new_action)

  def remove_branch(self):
    index, item = self.selection()
    if not item:
      self.statusBar().showMessage("No item selected.")
      return

    dialog = QtGui.QMessageBox.question(self, "Remove Dialog",
                "Remove %s and all children?" % item,
                buttons=QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok,
                defaultButton=QtGui.QMessageBox.Cancel)
    if dialog == QtGui.QMessageBox.Ok:
      for target in self.traverse(item):
        self.persistent_container.remove(target)

      if item.parent():
        item.parent().removeRow( item.row() )
      else:
        self.model.invisibleRootItem().removeRow( item.row() )

  def traverse(self, item, i=0):
    item_list = [item]
    while item.child(i):
      item_list.extend( self.traverse(item.child(i)) )
      i += 1
    return item_list

  def selection(self):
    if not self.view.selectedIndexes():
      return 0, None

    index = self.view.selectedIndexes()[0]
    item = self.model.itemFromIndex(index)

    return index, item

  def update(self):
    index, item = self.selection()
    self.settings.setWidget( item.module_widget if item else None )
    self.preview.setWidget( item.preview if item else None )

  def run_all(self):
    for item in self.traverse(self.model.invisibleRootItem())[1:]:
      item.module_widget.visible_widget.current_output()

  def save_project(self):
    print( "Save" )