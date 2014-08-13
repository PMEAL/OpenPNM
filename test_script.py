import numpy as np

import OpenPNM
import pytest

test_network = OpenPNM.Network.TestNet()

def test_print():
  print( test_network )

def test_add_boundary():
  pn = OpenPNM.Network.Cubic(name='net',loglevel=20,divisions=[5,5,5],
                             lattice_spacing=[0.0001],add_boundaries=True)

  keys_expected = {'pore.all', 'pore.back', 'pore.back_face', 'pore.bottom',
  'pore.bottom_face','pore.boundary', 'pore.coords', 'pore.front',
  'pore.front_face', 'pore.left', 'pore.left_face', 'pore.right',
  'pore.right_face', 'pore.top', 'pore.top_face', 'throat.all', 'throat.back',
  'throat.back_face', 'throat.bottom', 'throat.bottom_face', 'throat.boundary',
  'throat.conns', 'throat.front', 'throat.front_face', 'throat.left',
  'throat.left_face', 'throat.right', 'throat.right_face', 'throat.top',
  'throat.top_face'}

  keys_found = set(pn.keys())

  symmetric_diff = keys_found ^ keys_expected
  assert not symmetric_diff


def test_rectilinear_integrity():
  R = np.random.rand(2,3,4)
  # prune the mini way
  network = OpenPNM.Network.Template()
  network.generate(R)
  network.trim( ( R.ravel() <= R.mean() ) )
  O = network.asarray(network['pore.values'])
  # what it would look like normally
  M = np.where(R > R.mean(), R, 0)
  assert np.allclose(M, O)

if __name__ == '__main__':
  pytest.main()
