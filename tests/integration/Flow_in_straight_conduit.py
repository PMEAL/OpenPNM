# This test checks to ensure that pores which overlap give the same transport
# rates as pores which do not but are the same size
if __name__ == "__main__":

    import openpnm as op
    import numpy as np

    diameters = [0.9, 1.0, 1.1, 1.5]
    for D in diameters:
        pn5 = op.network.Cubic(shape=[5, 1, 1], spacing=1)
        pn5['pore.diameter'] = D
        pn5['throat.diameter'] = D
        pn5.regenerate_models()
        pn5.add_model(
            propname='throat.hydraulic_size_factors',
            model=op.models.geometry.hydraulic_size_factors.spheres_and_cylinders)

        pn3 = op.network.Cubic(shape=[3, 1, 1], spacing=2)
        pn3['pore.diameter'] = D
        pn3['throat.diameter'] = D
        pn3.add_model(
            propname='throat.hydraulic_size_factors',
            model=op.models.geometry.hydraulic_size_factors.spheres_and_cylinders)
        pn3.regenerate_models()

        ph1 = op.phase.Phase(network=pn3)
        ph1['pore.viscosity'] = 1.0
        ph1.add_model(
            propname='throat.hydraulic_conductance',
            model=op.models.physics.hydraulic_conductance.generic_hydraulic)

        ph2 = op.phase.Phase(network=pn5)
        ph2['pore.viscosity'] = 1.0
        ph2.add_model(
            propname='throat.hydraulic_conductance',
            model=op.models.physics.hydraulic_conductance.generic_hydraulic)

        sf1 = op.algorithms.StokesFlow(network=pn3, phase=ph1)
        sf1.set_value_BC(pores=pn3.pores('xmin'), values=1.0)
        sf1.set_value_BC(pores=pn3.pores('xmax'), values=0.0)
        sf1.run()
        r1 = sf1.rate(pores=pn3.pores('xmin'))

        sf2 = op.algorithms.StokesFlow(network=pn5, phase=ph1)
        sf2.set_value_BC(pores=pn5.pores('xmin'), values=1.0)
        sf2.set_value_BC(pores=pn5.pores('xmax'), values=0.0)
        sf2.run()
        r2 = sf2.rate(pores=pn5.pores('xmin'))
        assert np.allclose(r1, r2, atol=0, rtol=1e-10)
