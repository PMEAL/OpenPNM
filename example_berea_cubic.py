import openpnm as op
import numpy as np
import scipy as sp
import openpnm.models as mods
import matplotlib.pyplot as plt


def Pc_slit(target,
            throat_height='throat.height',
            throat_width='throat.width',
            contact_angle='throat.contact_angle',
            surface_tension='throat.surface_tension'):

    project = target.project
    network = project.network
    phase = project.find_phase(target)

    Ht = network[throat_height]
    Wt = network[throat_width]
    rt = 1/Ht + 1/Wt
    theta = phase[contact_angle]
    sigma = phase[surface_tension]
    Pc_slit = -sigma*sp.cos(sp.deg2rad(theta))*rt

    return Pc_slit[phase.throats(target.name)]


#def diffusive_conductance_slit(target,
#                               throat_height='throat.height',
#                               throat_width='throat.width',
#                               throat_length='throat.length',
#                               pore_height='pore.size_z',
#                               pore_width='pore.size_y',
#                               pore_length='pore.size_x',
#                               throat_diffusivity ='throat.diffusivity'):
#    project = target.project
#    net=project.network
#    conns=net['throat.conns']
#    phase=project.find_phase(target)
#    gd=sp.zeros((net.Nt, 3), dtype=float)
#
#    # Start with x-directional throats
#    Ts = net.throats('dir_x')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size, 1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = Ht.T*Wt
#    phase = water
#    gdt = (phase[throat_diffusivity][Ts]*Ta)/(Lt)
#    Hp = net['pore.size_z'][conns][Ts]
#    Lp = net['pore.size_x'][conns][Ts]/2
#    Wp = net['pore.size_y'][conns][Ts]
#    Pa = Hp*Wp
#    gdp1, gdp2 = ((phase[throat_diffusivity][conns][Ts]*Pa)/(Lp)).T
#    gd[Ts, :] = sp.vstack((gdp1, gdt, gdp2)).T
#    # y-directional throats
#    Ts = net.throats('dir_y')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size, 1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = Ht.T*Wt
#    gdt = (phase[throat_diffusivity][Ts]*Ta)/(Lt)
#    Hp = net['pore.size_z'][conns][Ts]
#    Lp = net['pore.size_y'][conns][Ts]/2
#    Wp = net['pore.size_x'][conns][Ts]
#    Pa = Hp*Wp
#    gdp1, gdp2 = ((phase[throat_diffusivity][conns][Ts]*Pa)/(Lp)).T
#    gd[Ts, :] = sp.vstack((gdp1, gdt, gdp2)).T
#    # z-directional throats
#    Ts = net.throats('dir_z')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size, 1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = Ht.T*Wt
#    gdt = (phase[throat_diffusivity][Ts]*Ta)/(Lt)
#    Hp = net['pore.size_x'][conns][Ts]
#    Lp = net['pore.size_z'][conns][Ts]/2
#    Wp = net['pore.size_y'][conns][Ts]
#    Pa = Hp*Wp
#    gdp1, gdp2 = ((phase[throat_diffusivity][conns][Ts]*Pa)/(Lp)).T
#    gd[Ts, :] = sp.vstack((gdp1, gdt, gdp2)).T
#    gdtotal=1/(sp.sum(1/gd, axis=1))
#    return gdtotal[phase.throat(target.name)]


def electric_conductance_slit(target, throat_height='throat.height',
                              throat_width='throat.width',
                              throat_length='throat.length',
                              pore_height='pore.size_z',
                              pore_width='pore.size_y',
                              pore_length='pore.size_x',
                              throat_electrical_conductivity='throat.electrical_conductivity'):

    project = target.project
    net = project.network
    conns = net['throat.conns']
    phase = project.find_phase(target)
    ge = sp.zeros((net.Nt, 3), dtype=float)
    phase['throat.electrical_conductivity'] = 1

    def getHtWt():
        global Ts
        Ts = net.throats('dir_x')
        global Ht
        Ht = sp.reshape(net[throat_height][Ts], (Ts.size,1))
        Wt = net[throat_width][Ts]
        Lt = net[throat_length][Ts]
        Ta = (Ht.T*Wt)
        global gte
        gte = (4.0*Ta*(phase[throat_electrical_conductivity][Ts]))/(Lt)
        return (Ts, Ht, Wt, Lt, Ta, gte)
    #throat_electrical_conductivity = np.full((15,15,15), 12.0,
    #                                      dtype =int)
    # Start with x-directional throats
    Ts = net.throats('dir_x')
    getHtWt()
    Hp = net['pore.size_z'][conns][Ts]
    #Lp = net['pore.size_x'][conns][Ts]/2
    Wp = net['pore.size_y'][conns][Ts]
    gpe1, gpe2 = ((2*np.pi*(Hp*Wp)*
                   phase[throat_electrical_conductivity][conns][Ts])/
                                    (1-np.log10(np.arcsin(Ht/Hp)))).T
    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
    # y-directional throats
    Ts = net.throats('dir_y')
    getHtWt()
    Hp = net['pore.size_z'][conns][Ts]
    #Lp = net['pore.size_y'][conns][Ts]/2
    Wp = net['pore.size_x'][conns][Ts]
    gpe1, gpe2 = ((2*np.pi*(Hp*Wp)*
                   phase[throat_electrical_conductivity][conns][Ts])/
                                   (1-np.log10(np.arcsin(Ht/Hp)))).T
    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
    # z-directional throats
    Ts = net.throats('dir_z')
    getHtWt()
    Hp = net['pore.size_x'][conns][Ts]
    #Lp = net['pore.size_z'][conns][Ts]/2
    Wp = net['pore.size_y'][conns][Ts]
    gpe1, gpe2 = ((2*np.pi*(Hp*Wp)*
                  phase[throat_electrical_conductivity][conns][Ts])/
                                    (1-np.log10(np.arcsin(Ht/Hp)))).T
    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
    getotal =1/(sp.sum(1/ge, axis=1))
    return getotal[phase.throats(target.name)]


#def electric_conductance_slit(target, throat_height='throat.height',
#                          throat_width='throat.width',
#                          throat_length='throat.length',
#                          pore_height='pore.size_z',
#                          pore_width='pore.size_y',
#                          pore_length='pore.size_x',
#                          throat_electrical_resistivity =
#                          'throat.electrical_resistivity'):
#
#    project = target.project
#    net=project.network
#    conns=net['throat.conns']
#    phase=project.find_phase(target)
#    ge=sp.zeros((net.Nt, 3), dtype=float)
#    phase['throat.electrical_resistivity'] =35.8
#    #throat_electrical_conductivity = np.full((15,15,15), 12.0,
#    #                                      dtype =int)
#    # Start with x-directional throats
#    Ts = net.throats('dir_x')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size,1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = (Ht.T*Wt)
#    gte = (4.0*Ta)/(Lt*(phase[throat_electrical_resistivity][Ts]))
#    Hp = net['pore.size_z'][conns][Ts]
#    #Lp = net['pore.size_x'][conns][Ts]/2
#    Wp = net['pore.size_y'][conns][Ts]
#    gpe1, gpe2 =((2*np.pi*(Hp*Wp))/(1-np.log10(np.arcsin(Ht/Hp))*
#                   phase[throat_electrical_resistivity][conns][Ts])).T
#    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
#    # y-directional throats
#    Ts = net.throats('dir_y')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size, 1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = (Ht.T*Wt)
#    gte = (4.0*Ta)/(Lt*(phase[throat_electrical_resistivity][Ts]))
#    Hp =net['pore.size_z'][conns][Ts]
#    #Lp = net['pore.size_y'][conns][Ts]/2
#    Wp = net['pore.size_x'][conns][Ts]
#    gpe1, gpe2 = ((2*np.pi*(Hp*Wp))/(1-np.log10(np.arcsin(Ht/Hp))*
#                   phase[throat_electrical_resistivity][conns][Ts])).T
#    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
#    # z-directional throats
#    Ts = net.throats('dir_z')
#    Ht = sp.reshape(net[throat_height][Ts], (Ts.size, 1))
#    Wt = net[throat_width][Ts]
#    Lt = net[throat_length][Ts]
#    Ta = (Ht.T*Wt)
#    gte = (4.0*Ta)/(Lt*(phase[throat_electrical_resistivity][Ts]))
#    Hp = net['pore.size_x'][conns][Ts]
#    #Lp = net['pore.size_z'][conns][Ts]/2
#    Wp = net['pore.size_y'][conns][Ts]
#    gpe1, gpe2 = ((2*np.pi*(Hp*Wp))/
#                                    (1-np.log10(np.arcsin(Ht/Hp))*
#                   phase[throat_electrical_resistivity][conns][Ts])).T
#    ge[Ts, :] = sp.vstack((gpe1, gte, gpe2)).T
#    getotal =1/(sp.sum(1/ge, axis=1))
#    return getotal[phase.throats(target.name)]


if __name__ == '__main__':
    proj = op.materials.BereaCubic(shape=[15, 15, 15])
    net = proj.network
    geo = proj[1]

    geo.add_model(propname='throat.size',
                  model=mods.geometry.throat_size.weibull,
                  shape=0.9, loc=1.1e-6, scale=0.000006,
                  seeds='throat.seed')
    geo.regenerate_models()
    plt.hist(x=geo['throat.size']*1e6, bins=25,
             weights=geo['throat.volume']*(1e3)**3, edgecolor='k')
    plt.hist(x=geo['pore.size_z']*1e6, bins=25,
             weights=geo['pore.volume']*(1e3)**3, edgecolor='k')

    # Below here, start to define phases and physics, then do simulations of
    # perm and pc curves for comparison to Marios paper.

    # Phase

    hg = op.phases.Mercury(network=net)
    air = op.phases.Air(network=net)
    hg['throat.contact_angle'] = 140
    hg['throat.surface_tension'] = 0.48

    phys_air = op.physics.GenericPhysics(network=net, phase=air,
                                         geometry=geo)
    phys_hg = op.physics.GenericPhysics(network=net, phase=hg,
                                        geometry=geo)

    phys_hg.add_model(propname='throat.entry_pressure', model=Pc_slit)

    mip = op.algorithms.Porosimetry(network=net, name='air')
    mip.setup(phase=hg)
    mip.set_inlets(net.pores(['top', 'bottom']))
    mip.run(points=25)
    mip.plot_intrusion_curve()
#    net.labels()
#    op.topotools.plot_connections(network =net)
#    op.topotools.plot_coordinates(network = net)

    water = op.phases.Water(network=net, name='water')
    water['throat.viscosity'] = 0.0001

    phys_water = op.physics.GenericPhysics(network=net, phase=water,
                                           geometry=geo)
    mod = op.models.physics.hydraulic_conductance.hagen_poiseuille_slit
    phys_water.add_model(propname='throat.hydraulic_conductance',
                         model=mod)
    phys_water.regenerate_models()

    alg = op.algorithms.StokesFlow(network=net, phase=water)
    BC1_pores = net.pores('pore.front')

    alg.set_value_BC(values=202650, pores=BC1_pores)

    BC2_pores = net.pores('pore.back')
    alg.set_value_BC(values=101325, pores=BC2_pores)
    alg.run()
    Q = alg.rate(pores=net.pores('front'))

    # Calculating porosity
    Vp = geo['pore.volume'][net.Ps]
    Vt = geo['throat.volume'][net.Ts]
    Vps = sp.sum(Vp)
    Vts = sp.sum(Vt)
    Vt = Vps + Vts
    Porosity = Vt/5.437e-9

    # Calculating permeability
    Lc = 0.0001
    A = (Lc**2)*15*15
    L = Lc*15
    mu = sp.mean(water['throat.viscosity'])

    Kxx = Q*mu*L/(A*101325)

    # Calculating electric current flow I
    water = op.phases.Water(network=net)
    #water['throat.viscosity'] = 0.001

    phys_water = op.physics.GenericPhysics(network=net, phase=water,
                                           geometry=geo)
    phys_water.add_model(propname='throat.electrical_conductance',
                         model=electric_conductance_slit)
    phys_water.regenerate_models()

    Om = op.algorithms.OhmicConduction(network=net, phase=water)
    BC1_pores = net.pores('pore.front')
    Om.set_value_BC(values=20, pores=BC1_pores)

    BC2_pores = net.pores('pore.back')
    Om.set_value_BC(values=0, pores=BC2_pores)
    Om.run()
    I = Om.rate(pores = net.pores('front'))

#    #Calculating Formation Factor F
#    Lc = 0.0001056
#    A = (Lc**2)*15*15
#    L = Lc*15
#    delV = 20
#    Rnet = delV*A/(I*L)
#    F = Rnet/1
#    # Calculating effective Diffusivity
#    phys_air.add_model(propname ='throat.diffusive_conductance',
#                       model=diffusive_conductance_slit)
#    fd = op.algorithms.FickianDiffusion(network=net, phase=air)
#    inlet = net.pores('front')
#    outlet = net.pores('back')
#    fd.set_value_BC(pores=inlet, values=0.8)
#    fd.set_value_BC(pores=outlet, values=0.0)
#    fd.run()
#    air.update(fd.results())
#    rate_inlet = fd.rate(pores=inlet)[0]
#    Kxx = alg.calc_effective_permeability(domain_area=A, domain_length=L)
#    water.update(alg.results())
#    alg.run()
#    # Calculating effective diffusivity
#    fd = op.algorithms.fickianDiffusion(network=net)
#    fd.setup(phase=air)
#    fd.set_
#    Na =fd.rate(pores = inlet)
#    C1 = 0.5
#    C2 = 0.02
#    Deff = Na*Lc/(A*(C1 - C2))
#
#    phys_air.add_model(model=diffusive_conductance_slit,
#                       propname='throat.diffusive_conductance')
#    phys_water.add_model(model=diffusive_conductance_slit,
#                         propname='throat.diffusive_conductance')
#    fd =op.algorithms.FickianDiffusion(network=net, phase =air)
#
#    fd.run()
#
#    phys_air.regenerate_models()
#    fd = op.algorithms.FickianDiffusion(network=net, phase=air)
#    fd.setup(phase=air)
#    inlet  = net.pores('front')
#    outlet = net.pores('back')
#    fd.set_value_BC(pores=inlet, values=1.0)
#    fd.set_value_BC(pores=outlet, values=0.0)
#    OP_1 = op.algorithms.OrdinaryPercolation(network=net)
#    OP_1.setup(phase=water, pore_volume='pore.volume',
#               throat_volume='throat.volume')
#    OP_1.set_inlets(net.pores('left'))
#    OP_1.run()
#
#    fig = OP_1.plot_intrusion_curve()
#    water.update(OP_1.results(Pc=10000))
#    air.update(OP_1.results(Pc=10000))
#
#
#    air.update(fd.results())
#    rate_inlet = fd.rate(pores=inlet)[0]