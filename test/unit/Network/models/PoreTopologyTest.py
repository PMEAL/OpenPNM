import OpenPNM
import scipy as sp


class PoreTopologyTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3], spacing=1)

    def test_get_subscripts(self):
        f = OpenPNM.Network.models.pore_topology.get_subscripts
        self.net.models.add(propname='pore.subscripts',
                            shape=[3, 3, 3],
                            model=f)
        assert 'pore.subscripts' in self.net.keys()
        assert sp.shape(self.net['pore.subscripts']) == (27, 3)

    def test_get_coords(self):
        f = OpenPNM.Network.models.pore_topology.adjust_spacing
        self.net.models.add(propname='pore.coords2',
                            model=f,
                            new_spacing=2)
        assert 'pore.coords2' in self.net.keys()
        a = sp.amax(self.net['pore.coords'])
        assert sp.amax(self.net['pore.coords2']) == 2*a
