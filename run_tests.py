import warnings
import OpenPNM
sim = OpenPNM.Base.Controller()
sim.loglevel = 50

def fxn():
    warnings.warn("future", FutureWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

import pytest
pytest.main()
