import scipy as sp
import matplotlib.pylab as plt

def Overview(net,
             throat_diameter='diameter',
             pore_diameter='diameter',
             throat_length='length',                         
             fig=None):
  r"""
  Plot a montage of key network size distribution histograms

  Parameters
  ----------
  net : OpenPNM Network Object
    The network for which the graphs are desired
  fig : Matplotlib figure object
    The canvas on which to draw the plots

  """
  if fig==None: fig = plt.figure()
  ax1 = fig.add_subplot(221)
  ax1.hist(net.get_pore_data(prop=pore_diameter)[net.get_pore_indices('internal')],25,facecolor='green')
  ax1.set_xlabel('Pore Diameter [m]')
  ax1.set_ylabel('Frequency')

  ax2 = fig.add_subplot(222)
  net.find_neighbor_pores(1)
  x = sp.zeros(net.num_pores())
  for i in list(range(0,sp.shape(net.adjacency_matrix['lil']['conns'].rows)[0])):
    x[i] = sp.shape(net.adjacency_matrix['lil']['conns'].rows[i])[0]
  ax2.hist(x,25,facecolor='yellow')
  ax2.set_xlabel('Coordination Number')
  ax2.set_ylabel('Frequency')

  ax3 = fig.add_subplot(223)
  ax3.hist(net.get_throat_data(prop=throat_diameter)[net.get_throat_indices('internal')],25,facecolor='blue')
  ax3.set_xlabel('Throat Diameter [m]')
  ax3.set_ylabel('Frequency')

  ax4 = fig.add_subplot(224)
  ax4.hist(net.get_throat_data(prop=throat_length)[net.get_throat_indices('internal')],25,facecolor='red')
  ax4.set_xlabel('Throat Length [m]')
  ax4.set_ylabel('Frequency')


def Capillary_Pressure_Curve(net,
                             fluid,
                             capillary_pressure='capillary_pressure',
                             pore_volume='volume',
                             throat_volume='volume',
                             fig=None):
  r"""
  Plot drainage capillary pressure curve

  Parameters
  ----------
  net : OpenPNM Network Object
      The network for which the graphs are desired
  fig : Matplotlib figure object
      Canvas on which to draw plots

  """
  if type(fluid)==str: fluid = net.find_object_by_name(fluid)
  try:
    PcPoints = sp.unique(fluid.get_throat_data(prop=capillary_pressure))
  except KeyError:
    raise Exception('Capillary pressure simulation has not been run')

  PcPoints = sp.unique(fluid.get_throat_data(prop=capillary_pressure))
  Snwp = sp.zeros_like(PcPoints)
  Ps = sp.r_[0:net.num_pores('internal')]
  for i in range(1,sp.size(PcPoints)):
      Pc = PcPoints[i]
      Snwp[i] = sum((fluid.get_throat_data(prop=capillary_pressure)[Ps]<Pc)*(net.get_throat_data(prop=throat_volume)[Ps]))/sum(net.get_throat_data(prop=throat_volume)[Ps])
  
  if fig==None: fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(PcPoints,Snwp,'r.-')
  ax.set_xlabel('Capillary Pressure')
  ax.set_ylabel('Fluid Saturation')
