import OpenPNM

def invading_fluid(
  name='water',
  Pc=2.206e6,
  Tc=647,
  MW=0.0181,
  ):
  diffusivity={'method': 'constant',
                'value': 1e-12}
  viscosity={'method': 'constant',
              'value': 0.001}
  molar_density={'method': 'constant',
                  'value': 44445}
  surface_tension={'method': 'Eotvos',
                        'k': 2.25e-4}
  contact_angle={'method': 'constant',
                  'value': 120}

  fluid = OpenPNM.Fluids.GenericFluid().create(locals())
  fluid.pore_conditions['temperature'] = 353
  fluid.pore_conditions['pressure'] = 101325
  return { 'invading_fluid': fluid }

def defending_fluid(
  name='air',
  Pc=3.771e6,
  Tc=132.65,
  MW=0.0291,
  ):
  diffusivity={'method': 'Fuller',
                   'MA': 0.03199,
                   'MB': 0.0291,
                   'vA': 16.3,
                   'vB': 19.7}
  viscosity={'method': 'Reynolds',
                 'uo': 0.001,
                  'b': 0.1}
  molar_density={'method': 'ideal_gas',
                      'R': 8.314}
  surface_tension={'method': 'constant',
                    'value': 0}
  contact_angle={'method': 'na'}

  fluid = OpenPNM.Fluids.GenericFluid().create(locals())
  fluid.pore_conditions['temperature'] = 353
  fluid.pore_conditions['pressure'] = 101325
  return { 'defending_fluid': fluid }




