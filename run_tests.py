import warnings

def fxn():
    warnings.warn("future", FutureWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

import pytest
pytest.main()
