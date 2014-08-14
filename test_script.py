import numpy as np

import OpenPNM
import pytest

def test_linear_solvers():
  pn = OpenPNM.Network.Cubic.empty(dims=[1,40,30], spacing=0.0001)
  geom = OpenPNM.Geometry.Toray090(network=pn,pores=pn.pores(),throats=pn.throats())
  air = OpenPNM.Phases.Air(network=pn)
  phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=pn.pores(),throats=pn.throats())

  BC1_pores = pn.pores(labels=['left'])
  BC2_pores = pn.pores(labels=['right'])

  alg_1 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_1.set_boundary_conditions(bctype='Dirichlet', bcvalue=1, pores=BC1_pores)
  alg_1.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_1.run(phase=air)

  alg_2 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_2.set_boundary_conditions(bctype='Neumann', bcvalue = -1e-11, pores=BC1_pores)
  alg_2.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_2.run(phase=air)

  alg_3 = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
  alg_3.set_boundary_conditions(bctype='Neumann_group', bcvalue = -3e-10, pores=BC1_pores)
  alg_3.set_boundary_conditions(bctype='Dirichlet', bcvalue=0, pores=BC2_pores)
  alg_3.run(phase=air)

  print( alg_1['pore.mole_fraction'][BC1_pores] )
  print( alg_2['pore.mole_fraction'][BC1_pores] )
  print( alg_3['pore.mole_fraction'][BC1_pores] )

def test_relative_pos_to_absolute_pos():
  pn = OpenPNM.Network.Cubic.empty([1,1,1], spacing=0.0001)
  np.testing.assert_array_equal( pn['pore.coords'], [0.5, 0.5, 0.5])

def incomplete_test_spacing_setting():
  pn = OpenPNM.Network.Cubic.empty([5,5,5], spacing=0.0001)
  print( pn['pore.coords'] )

def test_add_boundary():
  pn = OpenPNM.Network.Cubic.empty([5,5,5])
  pn.add_boundaries()

  keys_expected = {'pore.all', 'pore.back', 'pore.back_face', 'pore.bottom',
  'pore.bottom_face','pore.boundary', 'pore.coords', 'pore.front',
  'pore.front_face', 'pore.left', 'pore.left_face', 'pore.right',
  'pore.right_face', 'pore.top', 'pore.top_face', 'pore.values',
  'throat.all', 'throat.back', 'throat.back_face', 'throat.bottom',
  'throat.bottom_face', 'throat.boundary', 'throat.conns', 'throat.front',
  'throat.front_face', 'throat.left', 'throat.left_face', 'throat.right',
  'throat.right_face', 'throat.top', 'throat.top_face'}

  keys_found = set(pn.keys())

  symmetric_diff = keys_found ^ keys_expected
  assert not symmetric_diff

def test_rectilinear_integrity():
  R = np.random.rand(3,4,5)
  # prune the network way
  network = OpenPNM.Network.Cubic(R)
  network.trim( R.ravel() <= R.mean() )
  O = network.asarray(network['pore.values'])
  # what it would look like normally
  M = np.where(R > R.mean(), R, 0)
  assert np.allclose(M, O)

if __name__ == '__main__':
  pytest.main()
