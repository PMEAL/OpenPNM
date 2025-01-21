import numpy as np
import openpnm as op


d1 = {'pore.all': np.ones(10, dtype=bool),
      'throat.all': np.ones(10, dtype=bool),
      'param.test': 2.2,
      'phase1': {
          'pore.all': np.ones(10, dtype=bool),
          'pore.test1': np.ones(10, dtype=int),
          'throat.all': np.ones(20, dtype=bool),
          'throat.test1': np.ones(20, dtype=float),
          'param.test1': 2.2,
          'phase3': {
              'pore.all': np.ones(10, dtype=bool),
              'pore.test3': np.ones(10, dtype=int),
              'throat.all': np.ones(20, dtype=bool),
              'throat.test3': np.ones(20, dtype=float),
              'param.test3': 2.2,
              },
    },
      'phase2': {
          'pore.all': np.ones(10, dtype=bool),
          'pore.test2': np.ones(10, dtype=int),
          'throat.all': np.ones(20, dtype=bool),
          'throat.test2': np.ones(20, dtype=float),
          'param.test2': 2.2,
          'phase3': {
              'pore.all': np.ones(10, dtype=bool),
              'pore.test4': np.ones(10, dtype=int),
              'throat.all': np.ones(20, dtype=bool),
              'throat.test4': np.ones(20, dtype=float),
              'param.test4': 2.2,
              },
      },
}


d2 = op.pnmlib.core.flatten_dict(d1)
print('─'*10)
op.pnmlib.inspect.tree(op.pnmlib.core.get_data(target=d2, key='*.all'))
print('─'*10)
op.pnmlib.inspect.tree(op.pnmlib.core.get_data(target=d2, key='phase1/*.*'))
