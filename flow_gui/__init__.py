import sys
import logging
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(levelname)-8s :: %(message)-140s :: %(funcName)s[%(filename)s(%(lineno)d)]')

# ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
# ch.setFormatter(formatter)
# logger.addHandler(ch)

from PyQt4 import QtGui

from .pipeline import MainWindow

app = QtGui.QApplication(sys.argv)

main = MainWindow()
add_source = main.add_source
set_previews = main.set_previews

def run():
  main.show()
  app.exec_()
  app.deleteLater()