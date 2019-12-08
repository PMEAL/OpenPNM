import os
import openpnm as op
from pytest_notebook.nb_regression import NBRegressionFixture

fixture = NBRegressionFixture(exec_timeout=50)
fixture.diff_color_words = False

rootdir = os.path.split(os.path.split(op.__file__)[0])[0]
examples_dir = os.path.join(rootdir, r'examples')
test_dir = os.path.join(examples_dir, 'extractions')

for item in os.listdir(test_dir):
    if item.endswith('ipynb'):
        nbook = os.path.join(test_dir, item)
        fixture.check(str(nbook))
