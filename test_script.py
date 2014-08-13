import OpenPNM
import pytest

test_network = OpenPNM.Network.TestNet()

def test_print():
  print( test_network )

if __name__ == '__main__':
  pytest.main()
