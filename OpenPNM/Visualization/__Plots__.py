import scipy as sp

<<<<<<< HEAD
from __GenericVisualization__ import GenericVisualization

class Plots(GenericVisualization):
    r"""
    Methods for plotting common information about pore networks

    Parameters
    ----------

    Examples
    --------
    >>> print 'nothing yet'

    .. note::
    n/a

    """

    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(Plots,self).__init__(**kwargs)
        self._logger.debug("Construct Plots Class")

    @staticmethod
    def overview(net):
        r"""
        Plot a montage of key network size distribution histograms

        Parameters
        ----------
        net : OpenPNM Network Object
            The network for which the graphs are desired


        """
        import matplotlib.pyplot as plt
        #Print Pore and Throat size distributions
        plt.clf()
        plt.figure(1)
        plt.subplot(2,2,1)
        plt.hist(net.pore_properties['diameter'][net.pore_properties['type']==0],25,facecolor='green')
        plt.xlabel('Pore Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,2)
        net.get_neighbor_pores(1)
        x = sp.zeros(net.get_num_pores())
        for i in range(0,sp.shape(net.adjacency_matrix['lil']['connections'].rows)[0]):
            x[i] = sp.shape(net.adjacency_matrix['lil']['connections'].rows[i])[0]
        plt.hist(x,25,facecolor='yellow')
        plt.xlabel('Coordination Number')
        plt.ylabel('Frequency')
        plt.subplot(2,2,3)
        plt.hist(net.throat_properties['diameter'][net.throat_properties['type']==0],25,facecolor='blue')
        plt.xlabel('Throat Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,4)
        plt.hist(net.throat_properties['length'][net.throat_properties['type']==0],25,facecolor='red')
        plt.xlabel('Throat Length [m]')
        plt.ylabel('Frequency')

    @staticmethod
    def Capillary_Pressure_Curve(net):
        r"""
        Plot drainage capillary pressure curve

        Parameters
        ----------
        net : OpenPNM Network Object
            The network for which the graphs are desired

        """
        try:
            PcPoints = sp.unique(net.pore_conditions['Pc_invaded'])
        except:
            raise Exception('Capillary pressure simulation has not been run')
        import matplotlib.pyplot as plt
        PcPoints = sp.unique(net.pore_conditions['Pc_invaded'])
        Snwp = sp.zeros_like(PcPoints)
        Ps = sp.r_[0:net.get_num_pores([0])]
        for i in range(1,sp.size(PcPoints)):
            Pc = PcPoints[i]
            OpenPNM.Physics.MultiPhase.late_pore_filling(net,Pc=Pc)
            snwp = 1 - net.pore_conditions['satn_wp']
            Snwp[i] = sum(snwp[Ps]*(net.pore_conditions['Pc_invaded'][Ps]<Pc)*(net.pore_properties['volume'][Ps]))/sum(net.pore_properties['volume'][Ps])
        plt.plot(PcPoints,Snwp,'r.-')
        plt.xlabel('Capillary Pressure')
        plt.ylabel('Invading Fluid Saturation')
        plt.show(block = False)
=======
def Overview(net, fig=None):
  r"""
  Plot a montage of key network size distribution histograms

  Parameters
  ----------
  net : OpenPNM Network Object
    The network for which the graphs are desired
  fig : Matplotlib figure object
    The canvas on which to draw the plots

  """
  ax1 = fig.add_subplot(221)
  ax1.hist(net.pore_properties['diameter'][net.pore_properties['type']==0],25,facecolor='green')
  ax1.set_xlabel('Pore Diameter [m]')
  ax1.set_ylabel('Frequency')

  ax2 = fig.add_subplot(222)
  net.get_neighbor_pores(1)
  x = sp.zeros(net.get_num_pores())
  for i in range(0,sp.shape(net.adjacency_matrix['lil']['connections'].rows)[0]):
    x[i] = sp.shape(net.adjacency_matrix['lil']['connections'].rows[i])[0]
  ax2.hist(x,25,facecolor='yellow')
  ax2.set_xlabel('Coordination Number')
  ax2.set_ylabel('Frequency')

  ax3 = fig.add_subplot(223)
  ax3.hist(net.throat_properties['diameter'][net.throat_properties['type']==0],25,facecolor='blue')
  ax3.set_xlabel('Throat Diameter [m]')
  ax3.set_ylabel('Frequency')

  ax4 = fig.add_subplot(224)
  ax4.hist(net.throat_properties['length'][net.throat_properties['type']==0],25,facecolor='red')
  ax4.set_xlabel('Throat Length [m]')
  ax4.set_ylabel('Frequency')


def Capillary_Pressure_Curve(net, fig=None):
  r"""
  Plot drainage capillary pressure curve

  Parameters
  ----------
  net : OpenPNM Network Object
      The network for which the graphs are desired
  fig : Matplotlib figure object
      Canvas on which to draw plots

  """
  try:
    PcPoints = sp.unique(net.pore_conditions['Pc_invaded'])
  except KeyError:
    raise Exception('Capillary pressure simulation has not been run')
  
  PcPoints = sp.unique(net.pore_conditions['Pc_invaded'])
  Snwp = sp.zeros_like(PcPoints)
  Ps = sp.r_[0:net.get_num_pores([0])]
  for i in range(1,sp.size(PcPoints)):
      Pc = PcPoints[i]
      Snwp[i] = sum((net.pore_conditions['Pc_invaded'][Ps]<Pc)*(net.pore_properties['volume'][Ps]))/sum(net.pore_properties['volume'][Ps])
  
  ax = fig.add_subplot(111)
  ax.plot(PcPoints,Snwp,'r.-')
  ax.set_xlabel('Capillary Pressure')
  ax.set_ylabel('Invading Fluid Saturation')
>>>>>>> Weekend-Work
