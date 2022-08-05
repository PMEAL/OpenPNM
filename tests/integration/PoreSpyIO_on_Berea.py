if __name__ == "__main__":

    import openpnm as op
    import porespy as ps
    import numpy as np
    import os
    from pathlib import Path


    # %% Read image from file in fixtures
    path = Path(os.path.realpath(__file__),
                '../../../tests/fixtures/berea_100_to_300.npz')
    data = np.load(path.resolve())
    im = data['im']


    # %% Note meta data for this image
    data = {
        'shape': {
            'x': im.shape[0],
            'y': im.shape[1],
            'z': im.shape[2],
        },
        'resolution': 5.345e-6,
        'porosity': 19.6,
        'permeability': {
            'Kx': 1360,
            'Ky': 1304,
            'Kz': 1193,
            'Kave': 1286,
        },
        'formation factor': {
            'Fx': 23.12,
            'Fy': 23.99,
            'Fz': 25.22,
            'Fave': 24.08,
        },
    }

    # %% Perform extraction
    snow = ps.networks.snow2(im, voxel_size=data['resolution'],
                             boundary_width=[3, 0, 0], accuracy='standard')
    # ps.imshow(snow.regions/snow.phase)

    # %% Open network in OpenPNM
    pn = op.io.network_from_porespy(snow.network)
    pn['pore.diameter'] = pn['pore.equivalent_diameter']
    pn['throat.diameter'] = pn['throat.inscribed_diameter']
    pn['throat.spacing'] = pn['throat.total_length']
    # pn.add_model(propname='throat.conduit_lengths',
    #              model=op.models.geometry.conduit_lengths.pyramids_and_cuboids)
    pn.add_model(propname='throat.hydraulic_size_factors',
                 model=op.models.geometry.hydraulic_size_factors.pyramids_and_cuboids)
    pn.add_model(propname='throat.diffusive_size_factors',
                 model=op.models.geometry.diffusive_size_factors.pyramids_and_cuboids)

    pn.regenerate_models()
    h = op.utils.check_network_health(pn)
    op.topotools.trim(network=pn, pores=h['disconnected_pores'])
    gas = op.phase.Phase(network=pn)
    gas['pore.diffusivity'] = 1.0
    gas['pore.viscosity'] = 1.0
    gas.add_model_collection(op.models.collections.physics.basic)
    gas.regenerate_models()


    # %% Perform Fickian Diffusion to find formation factor
    fd = op.algorithms.FickianDiffusion(network=pn, phase=gas)
    fd.set_value_BC(pores=pn.pores('xmin'), values=1.0)
    fd.set_value_BC(pores=pn.pores('xmax'), values=0.0)
    fd.run()
    dC = 1.0
    L = (data['shape']['x'] + 6)*data['resolution']
    A = data['shape']['y']*data['shape']['z']*data['resolution']**2
    Deff = fd.rate(pores=pn.pores('xmin'))*(L/A)/dC
    F = 1/Deff
    print(f"The Formation factor of the extracted network is {F}")
    print(f"The compares to a value of {data['formation factor']['Fx']} from DNS")
    np.testing.assert_allclose(F, data['formation factor']['Fx'], rtol=0.09)

    # %% Perform Stokes flow to find Permeability coefficient
    sf = op.algorithms.StokesFlow(network=pn, phase=gas)
    sf.set_value_BC(pores=pn.pores('xmin'), values=1.0)
    sf.set_value_BC(pores=pn.pores('xmax'), values=0.0)
    sf.run()
    dP = 1.0
    L = (data['shape']['x'] + 6)*data['resolution']
    A = data['shape']['y']*data['shape']['z']*data['resolution']**2
    K = sf.rate(pores=pn.pores('xmin'))*(L/A)/dP*1e12
    print(f'Permeability coefficient is {K} Darcy')
    print(f"The compares to a value of {data['permeability']['Kx']/1000} from DNS")
    np.testing.assert_allclose(K, data['permeability']['Kx']/1000, rtol=0.05)
