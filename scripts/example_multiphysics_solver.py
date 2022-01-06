# %% Package imports
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
import openpnm.models.geometry.diffusive_size_factors as gd
np.random.seed(10)

# %% Set up for the solvers
Nx = 100
shape = [Nx, Nx, 1]
spacing = 1e-2  # 1/Nx adjust
net = op.network.Cubic(shape=shape, spacing=spacing)
# geo = op.geometry.CirclesAndRectangles(network=net, pores=net.Ps, throats=net.Ts)
# 2d square geometry
geo = op.geometry.GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)
geo.add_model(propname='pore.diameter', model=op.models.misc.constant, value=spacing)
geo.add_model(propname='throat.length', model=op.models.misc.constant, value=1e-15)
geo.add_model(propname='throat.diameter', model=op.models.misc.constant, value=spacing)
Ndim = op.topotools.dimensionality(network=net).sum()
Ap = At = spacing**(Ndim-1)
Vp = spacing**Ndim
geo["pore.cross_sectional_area"] = Ap
geo["throat.cross_sectional_area"] = At
geo["pore.volume"] = Vp
mod = gd.squares_and_rectangles if Ndim == 2 else gd.cubes_and_cuboids
geo.add_model(propname='throat.diffusive_size_factors', model=mod)

# Diff = op.metrics.EffectiveDiffusivity(network=net)
# # Diff.settings._update({
# # 'area': (10*1e-4)**(Ndim-1)})
# #'length':9*1e-4,})
# D_eff = Diff.run()
# print(f"Effective diffusivity in x direction is: {D_eff}")

air = op.phase.Air(network=net)
phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo)

# Make diffusivity a function of temperature
air['pore.temperature'] = 300

air.add_model(propname='pore.diffusivity',
              model=op.models.misc.linear,
              m=1.860793056e-06,
              b=-0.0005375624384,
              prop='pore.temperature')
phys.add_model(propname='throat.diffusive_conductance',
               model=op.models.physics.diffusive_conductance.generic_diffusive)

air.remove_model(propname='pore.thermal_conductivity')
air["pore.thermal_conductivity"] = 0.0262
phys.add_model(propname='throat.thermal_conductance',
               model=op.models.physics.thermal_conductance.generic_thermal)

tfd_settings = {
    "conductance": "throat.diffusive_conductance",
    "quantity": "pore.concentration",
}
tfc_settings = {
    "conductance": "throat.thermal_conductance",
    "quantity": "pore.temperature",
    "pore_volume": "pore.heat_capacity",
}

pardiso = op.solvers.PardisoSpsolve()
rk45 = op.integrators.ScipyRK45(verbose=True)

# %% Define Algorithms
# First algorithm, transient Fourier conduction
tfc = op.algorithms.TransientReactiveTransport(network=net, phase=air)
geo['pore.heat_capacity'] = geo['pore.volume'] * 1.0035 * 1000 * 1.225
tfc.settings._update(tfc_settings)
tfc.set_value_BC(net.pores("left"), 400)
T0 = np.ones(tfc.Np) * 300
T0[net.pores('left')] = 400

# Second algorithm, transient Fickian diffusion
tfd = op.algorithms.TransientReactiveTransport(network=net, phase=air)
# tfd.settings['variable_props'] = 'pore.temperature'
tfd.settings._update(tfd_settings)
tfd.set_value_BC(net.pores("left"), 100)
# tfd.set_value_BC(net.pores("right"), 100)
c0 = np.ones(tfd.Np)*50
c0[net.pores('left')] = 100

# Integrator parameters
t_initial = 0
t_final = 50
t_step = 5
n_steps = int((t_final - t_initial)/t_step) + 1
t = np.linspace(t_initial, t_final, n_steps)

# Add variable props to algs
tfd.set_variable_props('pore.temperature')
tfc.set_variable_props('pore.concentration')

# %% For Loop Solution
# Solve multiphysics system assuming temperature change is small over t_step
not_BC_pores = net.pores("left", mode='nor')
t_prev = 0
C_1_avg = [c0[not_BC_pores].mean()]
for ti in t:
    if ti == t_initial:
        continue
    print('time:', ti, "s")
    tspan = [t_prev, ti]
    t_prev = ti
    tout = ti
    # Solve for temperature first... add while loop
    # Temperature dependent thermal conductivity
    sol_1 = tfc.run(x0=T0, tspan=tspan, integrator=rk45, saveat=tout)
    air.regenerate_models()  # update diffusivuty
    phys.regenerate_models()  # update diffusive conductance
    sol_2 = tfd.run(x0=c0, tspan=tspan, integrator=rk45, saveat=tout)
    # Update initial coniditions
    T0 = sol_1[:, 1]
    c0 = sol_2[:, 1]
    C_1_avg.append(c0[not_BC_pores].mean())

# Note: dual coupling, time step needs to be small, hamed benchamrk solution
T_1 = sol_1[:, -1]
C_1 = sol_2[:, -1]


# %% Build RHS manually
def _build_rhs(algs):
    """Time-derivative of the variables being solved (T and c): dydt = RHS"""

    def ode_func(t, y):
        # initialize rhs
        rhs = []
        for i, alg in enumerate(algs):
            # get x from y, assume alg.Np is same for all algs
            x = y[i*tfc.Np:(i+1)*tfc.Np]
            # store x onto algorithm
            alg.x = x
            # build A and b
            alg._update_iterative_props()
            alg._build_A()
            alg._build_b()
            alg._apply_BCs()
            alg._apply_sources()
            A = alg.A.tocsc()
            b = alg.b
            V = geo[alg.settings['pore_volume']]
            rhs_alg = np.hstack(-A.dot(x) + b)/V
            rhs = np.hstack((rhs, rhs_alg))
        return rhs

    return ode_func


# Call solve_ivp and pass function that builds rhs as dydt
algs = [tfc, tfd]
rhs = _build_rhs(algs)
T0 = np.ones(tfc.Np) * 300
T0[net.pores('left')] = 400
c0 = np.ones(tfd.Np)*50
c0[net.pores('left')] = 100

y0 = np.hstack((T0, c0))  # ICs must include boundary condition
tspan = [0, t_final]
rtol = 1e-5
tmp = op.contrib.TransientMultiPhysics(algorithms=[tfc, tfd], network=net)
sol = tmp.run(y0, tspan, saveat=t)

# plot
# T_2 = sol.y[0:net.Np, -1]
# C_2 = sol.y[net.Np:, -1]
# C_2_avg = sol.y[net.Np:, 0:][not_BC_pores].mean(axis=0)
T_2 = sol.T.transpose()[0:net.Np, -1]
C_2 = sol.T.transpose()[net.Np:, -1]
C_2_avg = sol.T.transpose()[net.Np:, 0:][not_BC_pores].mean(axis=0)

# %% Pure Fickian Diffusion and Pure Fourier Conduction
air.remove_model(propname='pore.diffusivity')
air["pore.diffusivity"] = 2.067547840000001e-04
phys.regenerate_models()
tspan = [t_initial, t_final]
sol_tfd = tfd.run(x0=c0, tspan=tspan, integrator=rk45, saveat=t)
sol_tfc = tfc.run(x0=T0, tspan=tspan, integrator=rk45, saveat=t)

C_pure = sol_tfd.T[:, 0:][:, not_BC_pores].mean(axis=1)
T_pure = sol_tfc.T[:, 0:][:, not_BC_pores].mean(axis=1)

# %% Hardcode COMSOL data and make plots
# hardcode data
time = np.linspace(0, 2000, 401)
C_3_avg = np.array([50.09077301, 51.11933996, 51.58125199, 51.93146642,
                    52.23775725, 52.50851029, 52.73263011, 52.95421029,
                    53.16594848, 53.35453602, 53.54312356, 53.71682787,
                    53.88245681, 54.04808576, 54.18975466, 54.33138083,
                    54.47300700, 54.61463317, 54.75625934, 54.88513431,
                    55.00849572, 55.13185713, 55.25521854, 55.37857995,
                    55.49677443, 55.60715496, 55.71753550, 55.82791604,
                    55.93829658, 56.04723935, 56.14294177, 56.23864419,
                    56.33434661, 56.43004903, 56.52575145, 56.62145387,
                    56.71715629, 56.81285871, 56.90856113, 57.00426355,
                    57.09426635, 57.17851860, 57.26277085, 57.34702309,
                    57.43127534, 57.51552758, 57.59977983, 57.68403207,
                    57.76828432, 57.85253656, 57.93678881, 58.01358390,
                    58.08952824, 58.16547258, 58.24141692, 58.31736126,
                    58.39330560, 58.46924994, 58.54519428, 58.62113862,
                    58.69708296, 58.77114793, 58.84077335, 58.91039877,
                    58.98002419, 59.04964961, 59.11927504, 59.18890046,
                    59.25852588, 59.32815130, 59.39777672, 59.46740215,
                    59.53164782, 59.59355750, 59.65546718, 59.71737685,
                    59.77928653, 59.84119621, 59.90310589, 59.96501556,
                    60.02692524, 60.08883492, 60.15074460, 60.21265427,
                    60.27456395, 60.33647363, 60.39838330, 60.46029298,
                    60.52220266, 60.58411234, 60.64602201, 60.70793169,
                    60.76984137, 60.82856621, 60.88406630, 60.93956639,
                    60.99506649, 61.05056658, 61.10606667, 61.16156677,
                    61.21706686, 61.27256695, 61.32806705, 61.38356714,
                    61.43906723, 61.49456733, 61.55006742, 61.60556751,
                    61.66106761, 61.71656770, 61.77206779, 61.82756789,
                    61.88306798, 61.93856808, 61.99263524, 62.04330304,
                    62.09397084, 62.14463863, 62.19530643, 62.24597423,
                    62.29664203, 62.34730983, 62.39797762, 62.44864542,
                    62.49931322, 62.54998102, 62.60064882, 62.65131661,
                    62.70198441, 62.75265221, 62.80332001, 62.85398781,
                    62.90465560, 62.95532340, 63.00599120, 63.05629501,
                    63.10317825, 63.15006149, 63.19694473, 63.24382796,
                    63.29071120, 63.33759444, 63.38447768, 63.43136091,
                    63.47824415, 63.52512739, 63.57201063, 63.61889386,
                    63.66577710, 63.71266034, 63.75954358, 63.80642681,
                    63.85331005, 63.90019329, 63.94707652, 63.99395976,
                    64.04084300, 64.08667473, 64.13228346, 64.17778596,
                    64.22318220, 64.26847221, 64.31365597, 64.35873348,
                    64.40370475, 64.44856978, 64.49332856, 64.53798110,
                    64.58252740, 64.62696745, 64.67130126, 64.71552882,
                    64.75965014, 64.80366522, 64.84757405, 64.89137664,
                    64.93507298, 64.97866308, 65.02214694, 65.06552455,
                    65.10879592, 65.15196104, 65.19501992, 65.23797256,
                    65.28081895, 65.32355910, 65.36619300, 65.40872067,
                    65.45114208, 65.49345725, 65.53566618, 65.57776887,
                    65.61976531, 65.66165551, 65.70343946, 65.74511717,
                    65.78668863, 65.82810556, 65.86940812, 65.91060180,
                    65.95168661, 65.99266253, 66.03352957, 66.07428773,
                    66.11493701, 66.15547741, 66.19590893, 66.23623157,
                    66.27644533, 66.31655020, 66.35654620, 66.39643332,
                    66.43621155, 66.47588091, 66.51544138, 66.55489298,
                    66.59423569, 66.63346952, 66.67259447, 66.71161054,
                    66.75051773, 66.78931605, 66.82800547, 66.86658602,
                    66.90505769, 66.94342048, 66.98167439, 67.01981941,
                    67.05785556, 67.09578282, 67.13360121, 67.17131071,
                    67.20891134, 67.24640308, 67.28378594, 67.32105992,
                    67.35822503, 67.39562223, 67.43296786, 67.47022322,
                    67.50738832, 67.54446315, 67.58144772, 67.61834202,
                    67.65514606, 67.69185983, 67.72848333, 67.76501657,
                    67.80145955, 67.83781226, 67.87407470, 67.91024688,
                    67.94632879, 67.98232044, 68.01822182, 68.05403294,
                    68.08975379, 68.12538438, 68.16092470, 68.19637476,
                    68.23173455, 68.26700407, 68.30218333, 68.33727232,
                    68.37227105, 68.40717952, 68.44199771, 68.47672565,
                    68.51136331, 68.54591072, 68.58036785, 68.61473472,
                    68.64901133, 68.68319767, 68.71729375, 68.75129955,
                    68.78521510, 68.81962338, 68.85403937, 68.88839691,
                    68.92269602, 68.95693669, 68.99111893, 69.02524272,
                    69.05930808, 69.09331500, 69.12726348, 69.16115352,
                    69.19498513, 69.22875830, 69.26247303, 69.29612932,
                    69.32972717, 69.36326658, 69.39674756, 69.43017010,
                    69.46353420, 69.49683986, 69.53008709, 69.56327588,
                    69.59640622, 69.62947813, 69.66249161, 69.69544664,
                    69.72834324, 69.76118140, 69.79396112, 69.82668240,
                    69.85934524, 69.89194965, 69.92449562, 69.95698315,
                    69.98941224, 70.02178290, 70.05409511, 70.08634889,
                    70.11854423, 70.14853704, 70.17827746, 70.20801788,
                    70.23775829, 70.26749871, 70.29723913, 70.32697955,
                    70.35671997, 70.38646038, 70.41620080, 70.44594122,
                    70.47568164, 70.50542206, 70.53516248, 70.56490289,
                    70.59464331, 70.62438373, 70.65412415, 70.68386457,
                    70.71360499, 70.74334540, 70.77308582, 70.80282624,
                    70.83256666, 70.86230708, 70.89204750, 70.92178791,
                    70.95152833, 70.98126875, 71.01100917, 71.04074959,
                    71.07049000, 71.10023042, 71.12997084, 71.15971126,
                    71.18945168, 71.21919210, 71.24893251, 71.27867293,
                    71.30841335, 71.33708684, 71.36563625, 71.39418566,
                    71.42273507, 71.45128448, 71.47983389, 71.50838330,
                    71.53693271, 71.56548212, 71.59403153, 71.62258094,
                    71.65113035, 71.67967977, 71.70822918, 71.73677859,
                    71.76532800, 71.79387741, 71.82242682, 71.85097623,
                    71.87952564, 71.90807505, 71.93662446, 71.96517387,
                    71.99372328, 72.02227269, 72.05082210, 72.07937151,
                    72.10792092, 72.13647033, 72.16501974, 72.19356915,
                    72.22211856, 72.25066797, 72.27921739, 72.30776680,
                    72.33631621, 72.36388787, 72.39134584, 72.41880380,
                    72.44626176, 72.47371972, 72.50117769, 72.52863565,
                    72.55609361])
C_pure_comsol = np.array([50.09520053, 51.79746967, 52.54559098, 53.14200791,
                          53.64788429, 54.09629663, 54.49515365, 54.86004782,
                          55.19502501, 55.50859507, 55.80578977, 56.08665242,
                          56.35426276, 56.60994216, 56.85560370, 57.09641780,
                          57.33070876, 57.55760524, 57.77794304, 57.99294392,
                          58.20166629, 58.40411092, 58.60457597, 58.80071895,
                          58.99253985, 59.18003868, 59.36321544, 59.54207012,
                          59.71729637, 59.89137568, 60.06207532, 60.22939528,
                          60.39333556, 60.55389617, 60.71107709, 60.86636833,
                          61.02103287, 61.17330097, 61.32317262, 61.47064783,
                          61.61572660, 61.75840892, 61.90075658, 62.04263523,
                          62.18299896, 62.32184777, 62.45918166, 62.59500063,
                          62.72930468, 62.86209381, 62.99336802, 63.12312731,
                          63.25137168, 63.37810112, 63.50331565, 63.62712365,
                          63.75162001, 63.87491638, 63.99701274, 64.11790911,
                          64.23760549, 64.35610186, 64.47339824, 64.58949462,
                          64.70439100, 64.81808738, 64.93058377, 65.04188016,
                          65.15197655, 65.26203145, 65.37244722, 65.48201339,
                          65.59072998, 65.69859697, 65.80561438, 65.91178220,
                          66.01710043, 66.12156907, 66.22518812, 66.32795759,
                          66.42987746, 66.53094775, 66.63116844, 66.73194382,
                          66.83224788, 66.93192024, 67.03096092, 67.12936990,
                          67.22714719, 67.32429279, 67.42080669, 67.51668891,
                          67.61193943, 67.70655826, 67.80054540, 67.89390084,
                          67.98677774, 68.07943422, 68.17157436, 68.26320592,
                          68.35433665, 68.44497429, 68.53512659, 68.62480130,
                          68.71400616, 68.80274892, 68.89103734, 68.97887915,
                          69.06628211, 69.15325395, 69.23980244, 69.32593531,
                          69.41166032, 69.49698520, 69.58191772, 69.66646561,
                          69.75063662, 69.83443850, 69.91787900, 70.00096587,
                          70.08370684, 70.16610968, 70.24818212, 70.33009549,
                          70.41261687, 70.49483174, 70.57674009, 70.65834194,
                          70.73963726, 70.82062607, 70.90130837, 70.98168416,
                          71.06175343, 71.14151619, 71.22097243, 71.30012216,
                          71.37896537, 71.45750207, 71.53573226, 71.61365593,
                          71.69127309, 71.76858374, 71.84558787, 71.92228549,
                          71.99867659, 72.07476118, 72.15053926, 72.22601082,
                          72.30117587, 72.37603440, 72.45058642, 72.52545689,
                          72.60007095, 72.67442329, 72.74851390, 72.82234280,
                          72.89590997, 72.96921542, 73.04225915, 73.11504115,
                          73.18756144, 73.25982000, 73.33181684, 73.40355196,
                          73.47502536, 73.54623704, 73.61718699, 73.68787522,
                          73.75830173, 73.82846652, 73.89836958, 73.96801093,
                          74.03739055, 74.10650845, 74.17536463, 74.24395909,
                          74.31229182, 74.38036283, 74.44879865, 74.51714903,
                          74.58529153, 74.65322615, 74.72095288, 74.78847173,
                          74.85578270, 74.92288578, 74.98978098, 75.05646830,
                          75.12294774, 75.18921929, 75.25528296, 75.32113874,
                          75.38678665, 75.45222667, 75.51745880, 75.58248306,
                          75.64729943, 75.71190792, 75.77630852, 75.84050125,
                          75.90448609, 75.96826304, 76.03183212, 76.09519331,
                          76.15834662, 76.22129204, 76.28402958, 76.34655924,
                          76.40888102, 76.47099491, 76.53290092, 76.59459905,
                          76.65608929, 76.71737166, 76.77844613, 76.83931273,
                          76.89997144, 76.96042227, 77.02114300, 77.08177845,
                          77.14223416, 77.20251014, 77.26260638, 77.32252288,
                          77.38225965, 77.44181667, 77.50119397, 77.56039152,
                          77.61940933, 77.67824741, 77.73690575, 77.79538436,
                          77.85368322, 77.91180235, 77.96974174, 78.02750140,
                          78.08508131, 78.14248149, 78.19970194, 78.25674264,
                          78.31360361, 78.37028484, 78.42678633, 78.48310809,
                          78.53925011, 78.59521239, 78.65099493, 78.70659774,
                          78.76202081, 78.81726414, 78.87232773, 78.92721159,
                          78.98191571, 79.03644009, 79.09078474, 79.14494964,
                          79.19893482, 79.25274025, 79.30687487, 79.36096036,
                          79.41489608, 79.46868206, 79.52231827, 79.57580473,
                          79.62914143, 79.68232837, 79.73536556, 79.78825299,
                          79.84099066, 79.89357858, 79.94601674, 79.99830514,
                          80.05044379, 80.10243268, 80.15427181, 80.20596118,
                          80.25750080, 80.30889066, 80.36013077, 80.41122112,
                          80.46216171, 80.51295254, 80.56359362, 80.61408494,
                          80.66442650, 80.71461831, 80.76466036, 80.81455265,
                          80.86429518, 80.91388796, 80.96333099, 81.01262425,
                          81.06176776, 81.11076151, 81.15960550, 81.20829974,
                          81.25684422, 81.30523895, 81.35384545, 81.40239498,
                          81.45081606, 81.49910867, 81.54727282, 81.59530852,
                          81.64321576, 81.69099454, 81.73864485, 81.78616672,
                          81.83356012, 81.88082506, 81.92796155, 81.97496957,
                          82.02184914, 82.06860025, 82.11522290, 82.16171709,
                          82.20808282, 82.25432009, 82.30042891, 82.34640926,
                          82.39226116, 82.43798460, 82.48357958, 82.52904610,
                          82.57438416, 82.61959377, 82.66467491, 82.70962760,
                          82.75445183, 82.79914759, 82.84371490, 82.88815375,
                          82.93246415, 82.97664608, 83.02069956, 83.06462457,
                          83.10842113, 83.15208923, 83.19587681, 83.23959956,
                          83.28320845, 83.32670350, 83.37008469, 83.41335203,
                          83.45650551, 83.49954515, 83.54247093, 83.58528285,
                          83.62798093, 83.67056515, 83.71303552, 83.75539204,
                          83.79763470, 83.83976351, 83.88177847, 83.92367958,
                          83.96546683, 84.00714023, 84.04869978, 84.09014548,
                          84.13147732, 84.17269531, 84.21379944, 84.25478973,
                          84.29566616, 84.33642874, 84.37707747, 84.41761234,
                          84.45803336, 84.49834053, 84.53853384, 84.57861330,
                          84.61857891, 84.65843067, 84.69816858, 84.73779263,
                          84.77730283, 84.81669917, 84.85611216, 84.89544843,
                          84.93468278, 84.97381549, 85.01284686, 85.05177717,
                          85.09060670, 85.12933574, 85.16796457, 85.20649349,
                          85.24492276, 85.28325269, 85.32148356, 85.35961564,
                          85.39764923, 85.43558461, 85.47342207, 85.51116189,
                          85.54880435, 85.58634975, 85.62379837, 85.66115049,
                          85.69840640])
T_pure_comsol = np.array([300.17073109, 301.13246992, 301.61689048,
                          301.97339545, 302.28722624, 302.56111095,
                          302.79692630, 303.03274166, 303.24198241,
                          303.44291245, 303.63676187, 303.80006461,
                          303.96336735, 304.12667009, 304.28997283,
                          304.45231731, 304.59224339, 304.73216948,
                          304.87209556, 305.01202165, 305.15194773,
                          305.26895627, 305.38324198, 305.49752768,
                          305.61181338, 305.72609908, 305.84038479,
                          305.95467049, 306.06895619, 306.18324189,
                          306.29752759, 306.40216838, 306.50036407,
                          306.59855977, 306.69675547, 306.79495117,
                          306.89314687, 306.99134256, 307.08953826,
                          307.18773396, 307.28592966, 307.37869591,
                          307.45909393, 307.53949195, 307.61988997,
                          307.70028799, 307.78068601, 307.86108403,
                          307.94148205, 308.02188007, 308.10227809,
                          308.18267611, 308.26307413, 308.34347215,
                          308.42387017, 308.50426819, 308.58466621,
                          308.66506423, 308.74546225, 308.82586027,
                          308.90625829, 308.98665631, 309.05901376,
                          309.12818700, 309.19736024, 309.26653348,
                          309.33570672, 309.40487996, 309.47405320,
                          309.54322643, 309.61239967, 309.68157291,
                          309.75074615, 309.81991939, 309.88909263,
                          309.95826587, 310.02743911, 310.09661235,
                          310.16578559, 310.23495882, 310.30413206,
                          310.37330530, 310.44091679, 310.49784874,
                          310.55478069, 310.61171265, 310.66864460,
                          310.72557655, 310.78250850, 310.83944045,
                          310.89637241, 310.95330436, 311.01023631,
                          311.06716826, 311.12410021, 311.18103217,
                          311.23796412, 311.29489607, 311.35182802,
                          311.40875997, 311.46569193, 311.52262388,
                          311.57955583, 311.63648778, 311.69341973,
                          311.75035169, 311.80728364, 311.86421559,
                          311.92114754, 311.97807949, 312.03501145,
                          312.09194340, 312.14887535, 312.20580730,
                          312.26273925, 312.31967121, 312.37660316,
                          312.43353511, 312.49046706, 312.54739901,
                          312.60433097, 312.66126292, 312.71720030,
                          312.76633660, 312.81547291, 312.86460921,
                          312.91374552, 312.96288182, 313.01201813,
                          313.06115444, 313.11029074, 313.15942705,
                          313.20856335, 313.25769966, 313.30683596,
                          313.35597227, 313.40510857, 313.45424488,
                          313.50338119, 313.55251749, 313.60165380,
                          313.65079010, 313.69992641, 313.74906271,
                          313.79819902, 313.84733532, 313.89647163,
                          313.94560794, 313.99474424, 314.04388055,
                          314.09301685, 314.14215316, 314.19128946,
                          314.24042577, 314.28956207, 314.33869838,
                          314.38783468, 314.43697099, 314.48610730,
                          314.53524360, 314.58437991, 314.63351621,
                          314.68196102, 314.72567726, 314.76939349,
                          314.81310973, 314.85682596, 314.90054220,
                          314.94425844, 314.98797467, 315.03169091,
                          315.07540714, 315.11912338, 315.16283961,
                          315.20655585, 315.25027209, 315.29398832,
                          315.33770456, 315.38142079, 315.42513703,
                          315.46885326, 315.51256950, 315.55628574,
                          315.60000197, 315.64371821, 315.68743444,
                          315.73115068, 315.77486691, 315.81858315,
                          315.86229939, 315.90601562, 315.94973186,
                          315.99344809, 316.03716433, 316.08088056,
                          316.12459680, 316.16831304, 316.21202927,
                          316.25574551, 316.29946174, 316.34317798,
                          316.38689421, 316.43009868, 316.46980357,
                          316.50950846, 316.54921335, 316.58891824,
                          316.62862313, 316.66832802, 316.70803290,
                          316.74773779, 316.78744268, 316.82714757,
                          316.86685246, 316.90655735, 316.94626224,
                          316.98596713, 317.02567202, 317.06537691,
                          317.10508180, 317.14478669, 317.18449158,
                          317.22419647, 317.26390135, 317.30360624,
                          317.34331113, 317.38301602, 317.42272091,
                          317.46242580, 317.50213069, 317.54183558,
                          317.58154047, 317.62124536, 317.66095025,
                          317.70065514, 317.74036003, 317.78006491,
                          317.81976980, 317.85947469, 317.89917958,
                          317.93888447, 317.97858936, 318.01815040,
                          318.05669606, 318.09518553, 318.13361879,
                          318.17199586, 318.21031674, 318.24858141,
                          318.28678989, 318.32494217, 318.36303825,
                          318.40107813, 318.43906181, 318.47698930,
                          318.51486059, 318.55267568, 318.59043458,
                          318.62813727, 318.66578377, 318.70337407,
                          318.74090817, 318.77838608, 318.81580778,
                          318.85317329, 318.89048260, 318.92773571,
                          318.96493263, 319.00207335, 319.03915786,
                          319.07618619, 319.11315831, 319.15007424,
                          319.18693396, 319.22373749, 319.26048483,
                          319.29717596, 319.33381090, 319.37038963,
                          319.40691218, 319.44337852, 319.47978866,
                          319.51612384, 319.55227030, 319.58835324,
                          319.62437264, 319.66032851, 319.69622085,
                          319.73204966, 319.76781493, 319.80351668,
                          319.83915489, 319.87472957, 319.91024071,
                          319.94568833, 319.98107241, 320.01639297,
                          320.05164999, 320.08684348, 320.12197343,
                          320.15703986, 320.19204275, 320.22698211,
                          320.26185794, 320.29667024, 320.33141900,
                          320.36610423, 320.40072594, 320.43528411,
                          320.46977874, 320.50420985, 320.53857742,
                          320.57288147, 320.60712198, 320.64129896,
                          320.67541240, 320.70946232, 320.74344870,
                          320.77737155, 320.81123087, 320.84502666,
                          320.87875891, 320.91244968, 320.94623249,
                          320.97996039, 321.01363336, 321.04725141,
                          321.08081455, 321.11432276, 321.14777605,
                          321.18117442, 321.21451786, 321.24780639,
                          321.28104000, 321.31421868, 321.34734245,
                          321.38041129, 321.41342522, 321.44638422,
                          321.47928830, 321.51213746, 321.54493170,
                          321.57767102, 321.61035541, 321.64298489,
                          321.67555945, 321.70807908, 321.74054380,
                          321.77295359, 321.80530846, 321.83760841,
                          321.86985344, 321.90204355, 321.93417874,
                          321.96625901, 321.99828436, 322.03025478,
                          322.06217029, 322.09403087, 322.12583654,
                          322.15758728, 322.18928310, 322.22094224,
                          322.25267520, 322.28436036, 322.31599773,
                          322.34758730, 322.37912907, 322.41062305,
                          322.44206924, 322.47346763, 322.50481822,
                          322.53612102, 322.56737602, 322.59858323,
                          322.62974264, 322.66085425, 322.69191807,
                          322.72293410, 322.75390233, 322.78482276,
                          322.81569540, 322.84652024, 322.87729729,
                          322.90802654, 322.93870799, 322.96934165,
                          322.99992752, 323.03046559, 323.06095586,
                          323.09139834, 323.12179302, 323.15213991,
                          323.18243900, 323.21269029, 323.24289379,
                          323.27304950, 323.30315741, 323.33321752,
                          323.36322984, 323.39319436])

# Plot average concentrations comparing for loop, build rhs, and comsol
fig1, ax1 = plt.subplots()
# ax1.plot(t, C_1_avg, time, C_3_avg)
# ax1.legend(('For Loop', 'COMSOL'))
ax1.plot(t, C_1_avg, t, C_2_avg, time, C_3_avg)
ax1.legend(('For Loop', 'Build RHS', 'COMSOL'))
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Concentration mol/m3')
ax1.set_title('Average Concentration')

# plot average concentration from pure fickian diffusion
fig2, ax2 = plt.subplots()
ax2.plot(t, C_pure, time, C_pure_comsol)
ax2.legend(('OpenPNM', 'COMSOL'))
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Average Concentration (mol/m3)')
ax2.set_title('Only Diffusion')

# plot average concentration from pure fourier conduction
fig3, ax3 = plt.subplots()
ax3.plot(t, T_pure, time, T_pure_comsol)
ax3.legend(('OpenPNM', 'COMSOL'))
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Average Temperature (mol/m3)')
ax3.set_title('Only Conduction')

plt.figure(1)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_1.reshape((Nx, Nx)))
im_2 = ax[1].imshow(C_1.reshape((Nx, Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
im_1.set_clim(300, 400)
im_2.set_clim(0, 100)
fig.suptitle('For Loop', y=0.85)

plt.figure(2)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_2.reshape((Nx, Nx)))
im_2 = ax[1].imshow(C_2.reshape((Nx, Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
im_1.set_clim(300, 400)
im_2.set_clim(0, 100)
fig.suptitle('Build RHS', y=0.85)

# Error calculation
T_error = np.abs(T_2 - T_1) / T_2 * 100
C_error = np.abs(C_2 - C_1) / C_2 * 100
print(T_error.max(), C_error.max())

plt.figure(3)
fig, ax = plt.subplots(ncols=2)
im_1 = ax[0].imshow(T_error.reshape((Nx, Nx)))
im_2 = ax[1].imshow(C_error.reshape((Nx, Nx)))
fig.colorbar(im_1, ax=ax[0], fraction=0.046, pad=0.04)
fig.colorbar(im_2, ax=ax[1], fraction=0.046, pad=0.04)
ax[0].title.set_text('Temperature (K)')
ax[1].title.set_text('Concentration (mol m-3)')
plt.axis('off')
# im_1.set_clim(0, 1)
# im_2.set_clim(0, 1)
fig.suptitle('Error', y=0.85)

# %% Export Geometry
# op.io.COMSOL.export_data(network=net, filename='multiphysics_solver_2d')
# comparison for 2D square
