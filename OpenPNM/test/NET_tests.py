import unittest

class NETTest(unittest.TestCase):

  def setUp(self):
    import OpenPNM
    import scipy as sp

    self.sp = sp
    self.pn = OpenPNM.NET.BaseNet()

  def test_property_assignment(self):
    N = self.sp.rand(10,123)
    N2 = self.sp.rand(self.pn.get_num_pores(), 42)

    self.pn.pore_property(name='spam',  ndarray=N)
    self.pn.pore_property(name='monty', ndarray=N2)

    self.assertEqual(self.pn.get_num_pores(), 10)
    self.assertEqual(self.pn.get_num_throats(), 20)

    for key in self.pn.pore_properties.keys():
      self.assertEqual(self.pn.pore_properties[key].shape[0], 10)

    for key in self.pn.throat_properties.keys():
      self.assertEqual(self.pn.throat_properties[key].shape[0], 20)


if __name__ == '__main__':
  import __init__
  unittest.main()