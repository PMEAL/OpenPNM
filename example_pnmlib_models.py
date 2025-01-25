import openpnm.pnmlib as pnm
import numpy as np
import matplotlib.pyplot as plt
from openpnm.pnmlib.inspect import tree


# Create an empty project with just a plain dict
prj = {}

# Generate network topology
Nx, Ny, Nz = 2, 2, 1
pn = pnm.generators.cubic(
    shape=[Nx, Ny, Nz],
    node_prefix='pore',
    edge_prefix='throat',
    )

# Add labels to each half of the domain
pn['pore.left'] = pn['pore.coords'][:, 0] < 1
pn['pore.right'] = pn['pore.coords'][:, 0] > (Nx - 1)

# Not sure if this function is really necessary, but that's another story
pnm.core.add_network(prj, pn)

# Create a dictionary of models as follows:
# - The main key is the ultimate name of where the data will be stored once
#   generated. This can include domain notation (@) and group notation (/).
# - Each main key will hold a subdict which contains a handle to the model and
#   all the arguments the model needs.

onetime_models = {
    # We can put an array on network
    'network/pore.seed': {
        'model': pnm.models.geometry.random_seeds,
        'num_range': [0.2, 0.8],
    },
    # Or an array with the same name at the top level
    'pore.seed': {
        'model': pnm.models.geometry.random_seeds,
        'num_range': [0.2, 0.8],
    },
    # We can also assign scalar values to one domain
    'network/pore.volume@left': {
        'model': pnm.models.geometry.constant,
        'value': 3.3,
    },
    # And different values to another
    'network/pore.volume@right': {
        'model': pnm.models.geometry.constant,
        'value': 1.1,
    },
}

# We can have as many different model dicts as we want
geo_models1 = {
    # We can write to the top level using values from the network
    'pore.size3': {
        'model': pnm.models.geometry.product,
        'props': ['network/pore.seed', 'network/pore.seed'],
        },
    # We can use group and domain notation to write vectors.
    # Note that network/pore.seed is a full length array, so
    # the model also uses the @left when fetching values so we
    # don't need to include the @left in the props.
    'network/pore.size2@left': {
        'model': pnm.models.geometry.product,
        'props': ['network/pore.seed', 'network/pore.seed'],
    },
    'network/pore.size4@left': {
        'model': pnm.models.geometry.product,
        'props': ['network/pore.seed@left', 'network/pore.seed@left'],
    },
}

# We can "run" the models as follows
pnm.models.apply_models(prj, onetime_models)
pnm.models.apply_models(prj, geo_models1)

# This should also work with mixtures.  Let's define a phase with 2 components
pnm.core.set_data(prj, 'phase1/pore.temperature', 1.0)
pnm.core.set_data(prj, 'phase1/N2/param.MW', 28.0)
pnm.core.set_data(prj, 'phase1/N2/pore.x', 0.79)
pnm.core.set_data(prj, 'phase1/O2/param.MW', 32.0)
pnm.core.set_data(prj, 'phase1/O2/pore.x', 0.21)


# Here we have a model which finds the molecular weight of the mixture
def MW_mixture(target, MWs, xs):
    xs = [target[xs[i]] for i in range(len(xs))]
    MWs = [target[MWs[i]] for i in range(len(MWs))]
    MW_mix = np.sum([xs[i]*MWs[i] for i in range(len(xs))], axis=0)
    return MW_mix


mixture_models = {
    'phase1/pore.MW': {
        'model': MW_mixture,
        'MWs': ['phase1/O2/param.MW', 'phase1/N2/param.MW'],
        'xs': ['phase1/O2/pore.x', 'phase1/N2/pore.x'],
        },
    }
pnm.models.apply_models(prj, mixture_models)


# The above works, but is not generalizable. We would like to write the above
# 'model dict' without any reference to 'phase1' or 'O2' and 'N2'. We can avoid
# specifying 'phase1' in the model name by using the group name argument when
# calling apply_models, and rewriting the model to be totally generic.
def MW_mixture2(target, mixture, components=None, MW='param.MW', x='pore.x'):
    if components is None:
        components = pnm.core.get_components(target, mixture)
    MWs = [pnm.core.get_data(target, '/'.join((mixture, c, MW))) for c in components]
    xs = [pnm.core.get_data(target, '/'.join((mixture, c, x))) for c in components]
    MW_mix = np.sum([xs[i]*MWs[i] for i in range(len(xs))], axis=0)
    return MW_mix


# Instead of the mixture arg, which is more literal, we could use group, which is
# in keeping with the pnmlib terminology.  Technically, components are subgroups.
mixture_models = {
    'pore.MW2': {
        'model': MW_mixture2,
        'mixture': 'phase1',
        # 'components': ['N2', 'O2'],  # Looked up if not given
        'MW': 'param.MW',
        'x': 'pore.x',
        },
    }
# The group argument in the following function means we don't have to put the group
# in the model name above.  This allows the model to be generic.
pnm.models.apply_models(prj, mixture_models, group='phase1')


# I'm still not happy that we need to pass 'mixture' to the function. It would be
# truly generic if we could only get rid of that.
def MW_mixture3(target, components=None, MW='param.MW', x='pore.x'):
    if components is None:
        components = list({c.split('/')[0] for c in target.keys() if '/' in c})
    MWs = [pnm.core.get_data(target, '/'.join((c, MW))) for c in components]
    xs = [pnm.core.get_data(target, '/'.join((c, x))) for c in components]
    MW_mix = np.sum([xs[i]*MWs[i] for i in range(len(xs))], axis=0)
    return MW_mix


# Instead of the mixture arg, which is more literal, we could use group, which is
# in keeping with the pnmlib terminology.  Technically, components are subgroups.
mixture_models = {
    'pore.MW2': {
        'model': MW_mixture3,
        # 'components': ['N2', 'O2'],  # Looked up if not given
        'MW': 'param.MW',
        'x': 'pore.x',
        },
    }
# The group argument in the following function means we don't have to put the group
# in the model name above.  This allows the model to be generic.
pnm.models.apply_models(prj, mixture_models, group='phase1')

pnm.inspect.tree(prj)
pnm.inspect.info(prj)
