import pytest
import openpnm as op
import numpy as np


def random_seed(target, domain, seed=None, lim=[0, 1]):
    inds = target[domain]
    np.random.seed(seed)
    seeds = np.random.rand(inds.sum())*(lim[1]-lim[0]) + lim[0]
    return seeds


def factor(target, prop, f=1):
    vals = target[prop]*f
    return vals


def dolittle(target, domain):
    N = target[domain].sum()
    d = {}
    d['item1'] = np.ones([N, ])
    d['item2'] = np.ones([N, ])*2
    d['item3'] = np.ones([N, ])*3
    return d


class Base2Test:

    def setup_class(self):

        g = op.network.Cubic(shape=[3, 3, 1])

        g.add_model(propname='pore.seed@left',
                    model=random_seed,
                    lim=[0.2, 0.4],
                    regen_mode='deferred')
        g.add_model(propname='pore.seed@right',
                    model=random_seed,
                    lim=[0.7, 0.99],
                    regen_mode='deferred')
        g.add_model(propname='pore.seedx',
                    model=factor,
                    prop='pore.seed',
                    f=10,
                    regen_mode='deferred')
        g.add_model(propname='pore.test',
                    model=factor,
                    domain='front',
                    prop='pore.seed',
                    f=100,
                    regen_mode='deferred')
        g.add_model(propname='pore.dict2',
                    model=dolittle,
                    regen_mode='deferred',)
        g.add_model(propname='pore.dict3',
                    model=dolittle,
                    domain='left',
                    regen_mode='deferred',)
        g.add_model(propname='pore.dict3',
                    model=dolittle,
                    domain='right',
                    regen_mode='deferred',)

        # Use official args
        g.run_model('pore.seed', domain='pore.left')
        assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
        # Use partial syntax
        g.run_model('pore.seed', domain='left')
        assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
        # Use lazy syntax
        g.run_model('pore.seed@right')
        assert np.sum(np.isnan(g['pore.seed'])) == 3
        # Run the pore seed model for all domains at once:
        del g['pore.seed']
        g.run_model('pore.seed')  # This does not work yet
        assert 'pore.seed' in g.keys()
        assert g['pore.seed@left'].shape[0] == 3
        assert g['pore.seed@right'].shape[0] == 3
        assert 'pore.seedx' not in g.keys()
        # Full domain model
        g.run_model('pore.seedx')
        assert 'pore.seedx' in g.keys()
        x = g['pore.seedx']
        assert x[~np.isnan(x)].min() > 2
        # Run the non-domain enabled model on a subdomain
        assert 'pore.test' not in g
        g.run_model('pore.test')
        assert 'pore.test' in g
        assert np.isnan(g['pore.test']).sum() == 7
        del g['pore.test']
        g.run_model('pore.test@front')
        assert np.isnan(g['pore.test']).sum() == 7
        # Fetch data with lazy syntax
        assert g['pore.seed@left'].shape[0] == 3
        # Write data with lazy syntax, ensuring scalar to array conversion
        g['pore.seed@right'] = np.nan
        assert np.sum(~np.isnan(g['pore.seed'])) == g['pore.left'].sum()
        # Write array directly
        g['pore.seed@right'] = np.ones(3)*3
        assert np.sum(np.isnan(g['pore.seed'])) == 3
        # Use labels that were not used by models
        assert g['pore.seed@front'].shape[0] == 3
        # Write a dict
        g['pore.dict'] = {'item1': 1, 'item2': 2.0}
        assert g['pore.dict.item1'].sum() == 9
        assert g['pore.dict.item2'].sum() == 18
        # A dict with domains
        g['pore.dict'] = {'item1@left': 2, 'item2@right': 3.0}
        assert g['pore.dict.item1'].sum() == 12
        assert g['pore.dict.item2'].sum() == 21
        g['pore.dict'] = {'item3@left': 1, 'item3@right': 2.0}
        # Double dots...not sure how these should work
        g['pore.nested.name1'] = 10
        g['pore.nested.name2'] = 20
        assert isinstance(g['pore.nested'], dict)
        assert len(g['pore.nested']) == 2
        with pytest.raises(KeyError):
            g['pore.nested.fail']
        del g['pore.nested.name1']
        assert 'pore.nested.name1' not in g.keys()
        del g['pore.nested']
        assert 'pore.nested.name2' not in g.keys()
        # More fun
        c = g['conduit.seed']
        assert c.shape == (12, 3)
        assert 'throat.seed' not in g
        assert np.isnan(c[:, 1]).sum() == g.Nt
        # Run model that returns a dict to all pores
        g.run_model('pore.dict2')
        assert len(g['pore.dict2']) == 3
        assert g['pore.dict2.item1'].sum() == 9
        assert g['pore.dict2.item2'].sum() == 18
        # Run model that returns a dict to only subdomain pores
        g.run_model('pore.dict3@left')
        assert len(g['pore.dict3']) == 3
        assert np.isnan(g['pore.dict3.item1']).sum() == 6
        assert g['pore.dict3.item1@left'].sum() == 3
        assert np.isnan(g['pore.dict3.item2']).sum() == 6
        assert np.isnan(g['pore.dict3.item3']).sum() == 6
        # Run model on both domains that its assigned
        g.run_model('pore.dict3')
        assert g['pore.dict3.item1@left'].sum() == 3
        assert g['pore.dict3.item1@right'].sum() == 3


if __name__ == '__main__':

    t = Base2Test()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
