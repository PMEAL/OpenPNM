import os
import openpnm as op
import pytest
from pytest_notebook.nb_regression import NBRegressionFixture

# Initialize pytest-notebook 'fixture' class
fixture = NBRegressionFixture(exec_timeout=100)
fixture.diff_color_words = False

# Collect all notebooks in current directory
rootdir = os.path.split(os.path.split(op.__file__)[0])[0]
examples_dir = os.path.join(rootdir, 'examples')
test_dir = os.path.join(examples_dir, 'extractions')
skip = []


nbooks = []
for item in os.listdir(test_dir):
    if item.endswith('ipynb'):
        if item in skip:
            pass
        else:
            nbooks.append(os.path.join(test_dir, item))


@pytest.mark.parametrize("notebook", nbooks)
def test_notebook(notebook):
    fixture.check(str(notebook))


# def test_all_notebooks_in_folder(test_dir=test_dir, skip=skip):
#     for item in os.listdir(test_dir):
#         if item.endswith('ipynb'):
#             if item in skip:
#                 pass
#             else:
#                 nbook = os.path.join(test_dir, item)
#                 fixture.check(str(nbook))
