from PyQt4 import QtGui, QtCore

import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure


if __name__ == '__main__':

  app = QtGui.QApplication([])

  fig = Figure()
  ax = fig.add_axes([0,0,200,200])
  ax.plot( [1,2,3], [2,3,4] )
  
  f = FigureCanvas(fig)
  f.show()

  app.exec_()