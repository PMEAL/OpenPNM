import numpy as np

import OpenPNM
import pytest

test_network = OpenPNM.Network.TestNet()

def test_print():
  print( test_network )

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
