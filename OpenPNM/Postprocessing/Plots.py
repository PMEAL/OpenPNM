import scipy as _sp
import matplotlib.pylab as _plt

def distributions(net,
                 throat_diameter='throat.diameter',
                 pore_diameter='pore.diameter',
                 throat_length='throat.length'):
  r"""
  Plot a montage of key network size distribution histograms

  Parameters
  ----------
  net : OpenPNM Network Object
    The network for which the graphs are desired

  """
  fig = _plt.figure()
  ax1 = fig.add_subplot(221)
  ax1.hist(net[pore_diameter],25,facecolor='green')
  ax1.set_xlabel('Pore Diameter [m]')
  ax1.set_ylabel('Frequency')

  ax2 = fig.add_subplot(222)
  net.find_neighbor_pores(1)
  x = _sp.zeros(net.num_pores())
  for i in list(range(0,_sp.shape(net._adjacency_matrix['lil']['conns'].rows)[0])):
    x[i] = _sp.shape(net._adjacency_matrix['lil']['conns'].rows[i])[0]
  ax2.hist(x,25,facecolor='yellow')
  ax2.set_xlabel('Coordination Number')
  ax2.set_ylabel('Frequency')

  ax3 = fig.add_subplot(223)
  ax3.hist(net[throat_diameter],25,facecolor='blue')
  ax3.set_xlabel('Throat Diameter [m]')
  ax3.set_ylabel('Frequency')

  ax4 = fig.add_subplot(224)
  ax4.hist(net[throat_length],25,facecolor='red')
  ax4.set_xlabel('Throat Length [m]')
  ax4.set_ylabel('Frequency')


