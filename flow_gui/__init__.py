import sys
from PyQt4 import QtGui

from .pipeline import MainWindow

app = QtGui.QApplication(sys.argv)

main = MainWindow()
add_source = main.add_source

def run():
  main.show()
  app.exec_()
  app.deleteLater()