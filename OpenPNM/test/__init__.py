# allow tests to access files as if on main level
import sys
base_path = __file__.rsplit('/',3)[0]
# base_path = 
sys.path.append(base_path)
import OpenPNM

# test ALL suites
if __name__ == '__main__':
  import unittest

  test_modules = ['NET_tests']

  suite = unittest.TestSuite()

  for t in test_modules:
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))

  unittest.TextTestRunner().run(suite)