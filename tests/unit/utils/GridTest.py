import os
import pytest
import scipy as sp
import numpy as np
import openpnm as op
from pathlib import Path


class GridTest:

    def setup_class(self):
        self.ws = op.Workspace()
        self.ws.clear()
        self.proj = self.ws.new_project()
        self.net = op.network.Cubic(shape=[2, 2, 2], project=self.proj)
        Ps = self.net.pores('top')
        Ts = self.net.find_neighbor_throats(pores=Ps)
        self.geo1 = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                                throats=Ts)
        Ps = self.net.pores('bottom')
        Ts = ~self.net.to_mask(throats=Ts)
        self.geo2 = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                                throats=Ts)
        self.phase1 = op.phases.GenericPhase(network=self.net)
        self.phase2 = op.phases.GenericPhase(network=self.net)
        self.phys11 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase1,
                                                geometry=self.geo1)
        self.phys12 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase1,
                                                geometry=self.geo2)
        self.phys21 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase2,
                                                geometry=self.geo1)
        self.phys22 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase2,
                                                geometry=self.geo2)

    def test_add_and_remove_rows_and_cols(self):
        self.setup_class()
        g = self.proj.grid
        assert g.shape == (3, 3)
        g.add_col()
        assert g.shape == (3, 4)
        g.drop_col(3)
        assert g.shape == (3, 3)
        g.add_row()
        assert g.shape == (4, 3)
        g.drop_row(3)
        assert g.shape == (3, 3)

    def test_add_row_and_col_at_given_positions(self):
        self.setup_class()
        g = self.proj.grid
        g.add_row(1)
        assert g[1] == ['---', '---', '---']
        g.add_col(1)
        assert g[:][1] == ['---', '---', '---', '---']

    def test_nrows_ncols(self):
        self.setup_class()
        g = self.proj.grid
        assert g.nrows == 3
        g.add_row()
        assert g.nrows == 4
        assert g.ncols == 3
        g.add_col()
        assert g.ncols == 4

    def test_size_and_nnz(self):
        self.setup_class()
        g = self.proj.grid
        assert g.size == 9
        assert g.nnz == 9
        p = op.phases.GenericPhase(network=self.net)
        g = self.proj.grid
        assert g.size == 12
        assert g.nnz == 10

    def test_changing_blank(self):
        self.setup_class()
        p = op.phases.GenericPhase(network=self.net)
        g = self.proj.grid
        r = g.row(2)
        assert '---' in r
        g.blank = '==='
        r = g.row(2)
        assert '===' in r

    def test_set_col_and_row(self):
        self.setup_class()
        g = self.proj.grid
        g.set_col(0, None)
        assert None in g.col(0)
        g.set_row(0, None)
        assert None in g.row(0)

    def test_index_and_header(self):
        self.setup_class()
        g = self.proj.grid
        assert g.header._grid.table_data == [g.row(0)]
        assert g.index._grid.table_data == [[i] for i in g.col(0)]

    def test_phases_and_geometries(self):
        self.setup_class()
        g = self.proj.grid
        assert g.phases() == ['phase_01', 'phase_02']
        assert g.geometries() == ['geo_01', 'geo_02']

    def test_get_col_and_row_by_number_and_name(self):
        self.setup_class()
        g = self.proj.grid
        a = g.get_col(0)
        b = g.get_col('net_01')
        assert a._grid.table_data == b._grid.table_data
        a = g.get_row(1)
        b = g.get_row('geo_01')
        assert a._grid.table_data == b._grid.table_data

    def test_set_row_and_col(self):
        self.setup_class()
        g = self.proj.grid
        g.set_row_and_col(row='x', col='z', val='z')
        assert g[-1][-1] == 'z'
        g.set_row_and_col(row='x', col='z', val='zz')
        assert g[-1][-1] == 'zz'

    def test_style_change(self):
        self.setup_class()
        g = self.proj.grid
        assert g.style == 'ascii'
        g.style = 'double'
        assert g.style == 'double'


if __name__ == '__main__':

    t = GridTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
