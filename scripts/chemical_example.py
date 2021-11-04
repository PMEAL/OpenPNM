import openpnm as op
from openpnm.phases import mixtures
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)

o2 = mixtures.GasByName(network=pn, species='oxygen', name='o2')
n2 = mixtures.GasByName(network=pn, species='nitrogen', name='n2')
air = mixtures.GasMixture(network=pn, components=[o2, n2])
air.set_mole_fraction(component=o2, values=0.21)
air.update_mole_fractions()
air.regenerate_models()

water = mixtures.LiquidByName(network=pn, species='water', name='h2o')
ethanol = mixtures.LiquidByName(network=pn, species='ethanol', name='etOH')
vodka = mixtures.LiquidMixture(network=pn, components=[water, ethanol])
vodka.set_mole_fraction(component=ethanol, values=0.40)
vodka.update_mole_fractions()
vodka.regenerate_models()
