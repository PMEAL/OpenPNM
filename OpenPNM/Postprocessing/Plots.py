import scipy as _sp
import matplotlib.pylab as _plt


def profiles(network, fig=None, values=None, bins=[10, 10, 10]):
    r"""
    Compute the profiles for the property of interest and plots it in all
    three dimensions

    Parameters
    ----------
    network : OpenPNM Network object

    values : array_like, optional
        The pore property values to be plotted as a profile

    bins : int or list of ints, optional
        The number of bins to divide the domain into for averaging.

    Notes
    -----
    Either propname or values can be sent, but not both

    """
    if fig is None:
        fig = _plt.figure()
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    ax = [ax1, ax2, ax3]
    xlab = ['x coordinate', 'y_coordinate', 'z_coordinate']
    for n in [0, 1, 2]:
        n_min, n_max = [_sp.amin(network['pore.coords'][:, n]),
                        _sp.amax(network['pore.coords'][:, n])]
        steps = _sp.linspace(n_min, n_max, bins[n]+1, endpoint=True)
        vals = _sp.zeros_like(steps)
        for i in range(0, len(steps)-1):
            temp = (network['pore.coords'][:, n] > steps[i]) * \
                   (network['pore.coords'][:, n] <= steps[i+1])
            vals[i] = _sp.mean(values[temp])
        yaxis = vals[:-1]
        xaxis = (steps[:-1] + (steps[1]-steps[0])/2)/n_max
        ax[n].plot(xaxis, yaxis, 'bo-')
        ax[n].set_xlabel(xlab[n])
        ax[n].set_ylabel('Slice Value')
    fig.show()


def porosity_profile(network, fig=None, axis=2):

    r"""
    Compute and plot the porosity profile in all three dimensions

    Parameters
    ----------
    network : OpenPNM Network object
    axis : integer type 0 for x-axis, 1 for y-axis, 2 for z-axis

    Notes
    -----
    the area of the porous medium at any position is calculated from the
    maximum pore coordinates in each direction

    """
    if fig is None:
        fig = _plt.figure()
    L_x = _sp.amax(network['pore.coords'][:, 0]) + \
        _sp.mean(((21/88.0)*network['pore.volume'])**(1/3.0))
    L_y = _sp.amax(network['pore.coords'][:, 1]) + \
        _sp.mean(((21/88.0)*network['pore.volume'])**(1/3.0))
    L_z = _sp.amax(network['pore.coords'][:, 2]) + \
        _sp.mean(((21/88.0)*network['pore.volume'])**(1/3.0))
    if axis is 0:
        xlab = 'x-direction'
        area = L_y*L_z
    elif axis is 1:
        xlab = 'y-direction'
        area = L_x*L_z
    else:
        axis = 2
        xlab = 'z-direction'
        area = L_x*L_y
    n_max = _sp.amax(network['pore.coords'][:, axis]) + \
        _sp.mean(((21/88.0)*network['pore.volume'])**(1/3.0))
    steps = _sp.linspace(0, n_max, 100, endpoint=True)
    vals = _sp.zeros_like(steps)
    p_area = _sp.zeros_like(steps)
    t_area = _sp.zeros_like(steps)

    rp = ((21/88.0)*network['pore.volume'])**(1/3.0)
    p_upper = network['pore.coords'][:, axis] + rp
    p_lower = network['pore.coords'][:, axis] - rp
    TC1 = network['throat.conns'][:, 0]
    TC2 = network['throat.conns'][:, 1]
    t_upper = network['pore.coords'][:, axis][TC1]
    t_lower = network['pore.coords'][:, axis][TC2]

    for i in range(0, len(steps)):
        p_temp = (p_upper > steps[i])*(p_lower < steps[i])
        t_temp = (t_upper > steps[i])*(t_lower < steps[i])
        p_area[i] = sum((22/7.0)*(rp[p_temp]**2 -
                        (network['pore.coords'][:, axis][p_temp]-steps[i])**2))
        t_area[i] = sum(network['throat.area'][t_temp])
        vals[i] = (p_area[i]+t_area[i])/area
    yaxis = vals
    xaxis = steps/n_max
    _plt.plot(xaxis, yaxis, 'bo-')
    _plt.xlabel(xlab)
    _plt.ylabel('Porosity')
    fig.show()


def saturation_profile(network, phase, fig=None, axis=2):

    r"""
    Compute and plot the saturation profile in all three dimensions

    Parameters
    ----------
    network : OpenPNM Network object
    phase : the invading or defending phase to plot its saturation distribution
    axis : integer type 0 for x-axis, 1 for y-axis, 2 for z-axis

    """
    if fig is None:
        fig = _plt.figure()
    if phase is None:
        raise Exception('The phase for saturation profile plot is not given')
    if axis is 0:
        xlab = 'x-direction'
    elif axis is 1:
        xlab = 'y-direction'
    else:
        axis = 2
        xlab = 'z-direction'
    n_max = _sp.amax(network['pore.coords'][:, axis]) + \
        _sp.mean(((21/88.0)*network['pore.volume'])**(1/3.0))
    steps = _sp.linspace(0, n_max, 100, endpoint=True)
    p_area = _sp.zeros_like(steps)
    op_area = _sp.zeros_like(steps)
    t_area = _sp.zeros_like(steps)
    ot_area = _sp.zeros_like(steps)
    vals = _sp.zeros_like(steps)
    PO = phase['pore.occupancy']
    TO = phase['throat.occupancy']
    rp = ((21/88.0)*network['pore.volume'])**(1/3.0)
    p_upper = network['pore.coords'][:, axis] + rp
    p_lower = network['pore.coords'][:, axis] - rp

    TC1 = network['throat.conns'][:, 0]
    TC2 = network['throat.conns'][:, 1]
    t_upper = network['pore.coords'][:, axis][TC1]
    t_lower = network['pore.coords'][:, axis][TC2]

    for i in range(0, len(steps)):
        op_temp = (p_upper > steps[i])*(p_lower < steps[i])*PO
        ot_temp = (t_upper > steps[i])*(t_lower < steps[i])*TO
        op_temp = _sp.array(op_temp, dtype='bool')
        ot_temp = _sp.array(op_temp, dtype='bool')
        p_temp = (p_upper > steps[i])*(p_lower < steps[i])
        t_temp = (t_upper > steps[i])*(t_lower < steps[i])
        op_area[i] = sum((22/7.0)*(rp[op_temp]**2 -
                         (network['pore.coords'][:, axis][op_temp]-steps[i])**2))
        ot_area[i] = sum(network['throat.area'][ot_temp])
        p_area[i] = sum((22/7.0)*(rp[p_temp]**2 -
                        (network['pore.coords'][:, axis][p_temp]-steps[i])**2))
        t_area[i] = sum(network['throat.area'][t_temp])
        vals[i] = (op_area[i]+ot_area[i])/(p_area[i]+t_area[i])
        if vals[i] > 1:
            vals[i] = 1.
        if _sp.isnan(vals[i]):
            vals[i] = 1.

    if vals[-1] == 1.:
        vals = vals[::-1]

    yaxis = vals
    xaxis = steps/n_max
    _plt.plot(xaxis, yaxis, 'bo-')
    _plt.xlabel(xlab)
    _plt.ylabel('Saturation')
    fig.show()


def distributions(net, fig=None, throat_diameter='throat.diameter',
                  pore_diameter='pore.diameter', throat_length='throat.length',
                  exclude_boundaries=True, geom_list=None):
    r"""
    Plot a montage of key network size distribution histograms

    Parameters
    ----------
    net : OpenPNM Network Object
    The network for which the graphs are desired

    """
    if fig is None:
        fig = _plt.figure()

    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(wspace=0.4)

    if geom_list is not None:
        include_pores = [False]*net.num_pores()
        include_throats = [False]*net.num_throats()
        for geom in geom_list:
            include_pores = include_pores | net['pore.' + geom]
            include_throats = include_throats | net['throat.' + geom]
    else:
        include_pores = net['pore.all']
        include_throats = net['throat.all']
    pores = net.pores()[include_pores]
    throats = net.throats()[include_throats]

    ax1 = fig.add_subplot(221)
    ax1.hist(net[pore_diameter][pores], 25, facecolor='green')
    ax1.set_xlabel('Pore Diameter')
    ax1.set_ylabel('Frequency')
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax2 = fig.add_subplot(222)
    x = net.num_neighbors(pores, flatten=False)
    ax2.hist(x, 25, facecolor='yellow')
    ax2.set_xlabel('Coordination Number')
    ax2.set_ylabel('Frequency')

    ax3 = fig.add_subplot(223)
    ax3.hist(net[throat_diameter][throats], 25, facecolor='blue')
    ax3.set_xlabel('Throat Diameter')
    ax3.set_ylabel('Frequency')
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax4 = fig.add_subplot(224)
    ax4.hist(net[throat_length][throats], 25, facecolor='red')
    ax4.set_xlabel('Throat Length')
    ax4.set_ylabel('Frequency')
    ax4.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    fig.show()


def pore_size_distribution(network, fig=None):
    r"""
    Plot the pore and throat size distribution which is the accumulated
    volume vs. the diameter in a semilog plot

    Parameters
    ----------
    network : OpenPNM Network object

    """
    if fig is None:
        fig = _plt.figure()
    dp = network['pore.diameter']
    Vp = network['pore.volume']
    dt = network['throat.diameter']
    Vt = network['throat.volume']
    dmax = max(max(dp), max(dt))
    steps = _sp.linspace(0, dmax, 100, endpoint=True)
    vals = _sp.zeros_like(steps)
    for i in range(0, len(steps)-1):
        temp1 = dp > steps[i]
        temp2 = dt > steps[i]
        vals[i] = sum(Vp[temp1]) + sum(Vt[temp2])
    yaxis = vals
    xaxis = steps
    _plt.semilogx(xaxis, yaxis, 'b.-')
    _plt.xlabel('Pore & Throat Diameter (m)')
    _plt.ylabel('Cumulative Volume (m^3)')
    fig.show()


def drainage_curves(inv_alg, fig=None, Pc='inv_Pc', sat='inv_sat',
                    seq='inv_seq', timing=None):
    r"""
    Plot a montage of key saturation plots

    Parameters
    ----------
    inv_alg : OpenPNM Algorithm Object
    The invasion algorithm for which the graphs are desired

    timing : string
    if algorithm keeps track of simulated time, insert string here

    """
    inv_throats = inv_alg.toindices(inv_alg['throat.' + seq] > 0)
    sort_seq = _sp.argsort(inv_alg['throat.'+seq][inv_throats])
    inv_throats = inv_throats[sort_seq]

    if fig is None:
        fig = _plt.figure(num=1, figsize=(13, 10), dpi=80,
                          facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(231)   # left
    ax2 = fig.add_subplot(232)   # middle
    ax3 = fig.add_subplot(233)   # right
    ax4 = fig.add_subplot(234)   # left
    ax5 = fig.add_subplot(235)   # middle
    ax6 = fig.add_subplot(236)   # right

    ax1.plot(inv_alg['throat.' + Pc][inv_throats],
             inv_alg['throat.' + sat][inv_throats])
    ax1.set_xlabel('Capillary Pressure (Pa)')
    ax1.set_ylabel('Saturation')
    ax1.set_ylim([0, 1])
    ax1.set_xlim([0.99*min(inv_alg['throat.' + Pc][inv_throats]),
                  1.01*max(inv_alg['throat.' + Pc][inv_throats])])

    ax2.plot(inv_alg['throat.' + seq][inv_throats],
             inv_alg['throat.' + sat][inv_throats])
    ax2.set_xlabel('Simulation Step')
    ax2.set_ylabel('Saturation')
    ax2.set_ylim([0, 1])
    ax2.set_xlim([0, 1.01*max(inv_alg['throat.' + seq][inv_throats])])

    if timing is None:
        ax3.plot(0, 0)
        ax3.set_xlabel('No Time Data Available')
    else:
        ax3.plot(inv_alg['throat.' + timing][inv_throats],
                 inv_alg['throat.' + sat][inv_throats])
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Saturation')
        ax3.set_ylim([0, 1])
        ax3.set_xlim([0, 1.01*max(inv_alg['throat.' + timing][inv_throats])])

    ax4.plot(inv_alg['throat.' + sat][inv_throats],
             inv_alg['throat.' + Pc][inv_throats])
    ax4.set_ylabel('Capillary Pressure (Pa)')
    ax4.set_xlabel('Saturation')
    ax4.set_xlim([0, 1])
    ax4.set_ylim([0.99*min(inv_alg['throat.' + Pc][inv_throats]),
                  1.01*max(inv_alg['throat.' + Pc][inv_throats])])

    ax5.plot(inv_alg['throat.' + seq][inv_throats],
             inv_alg['throat.' + Pc][inv_throats])
    ax5.set_xlabel('Simulation Step')
    ax5.set_ylabel('Capillary Pressure (Pa)')
    ax5.set_ylim([0.99*min(inv_alg['throat.' + Pc][inv_throats]),
                  1.01*max(inv_alg['throat.' + Pc][inv_throats])])
    ax5.set_xlim([0, 1.01*max(inv_alg['throat.' + seq][inv_throats])])

    if timing is None:
        ax6.plot(0, 0)
        ax6.set_xlabel('No Time Data Available')
    else:
        ax6.plot(inv_alg['throat.' + timing][inv_throats],
                 inv_alg['throat.' + Pc][inv_throats])
        ax6.set_xlabel('Time (s)')
        ax6.set_ylabel('Capillary Pressure (Pa)')
        ax6.set_ylim([0.99*min(inv_alg['throat.' + Pc][inv_throats]),
                      1.01*max(inv_alg['throat.' + Pc][inv_throats])])
        ax6.set_xlim([0, 1.01*max(inv_alg['throat.' + timing][inv_throats])])

    fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.1)
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)
    ax5.grid(True)
    ax6.grid(True)
    fig.show()
