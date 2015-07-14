import OpenPNM
import OpenPNM.Utilities.vertexops as vo


class DelaunayCubicTest:

    def test_simple_cubic(self):
        pn = OpenPNM.Network.DelaunayCubic(shape=[3, 3, 3])
        pn.add_boundaries()
        B1 = pn.pores("top_boundary")
        B2 = pn.pores("bottom_boundary")
        assert vo.vertex_dimension(pn, B1, B2, 'minmax') == [0.0, 3.0, 0.0,
                                                             3.0, 0.0, 3.0]
        assert pn.num_pores() == 81

    def test_orthorhombic(self):
        pn = OpenPNM.Network.DelaunayCubic(shape=[3, 3, 3], arrangement='O')
        pn.add_boundaries()
        B1 = pn.pores("top_boundary")
        B2 = pn.pores("bottom_boundary")
        assert vo.vertex_dimension(pn, B1, B2, 'minmax') == [0.0, 3.0, 0.0,
                                                             3.0, 0.0, 3.0]
        assert pn.num_pores() == 81

    def test_bcc(self):
        pn = OpenPNM.Network.DelaunayCubic(shape=[3, 3, 3], arrangement='BCC')
        pn.add_boundaries()
        B1 = pn.pores("top_boundary")
        B2 = pn.pores("bottom_boundary")
        assert vo.vertex_dimension(pn, B1, B2, 'minmax') == [0.0, 3.0, 0.0,
                                                             3.0, 0.0, 3.0]
        assert pn.num_pores() == 89

    def test_fcc(self):
        pn = OpenPNM.Network.DelaunayCubic(shape=[3, 3, 3], arrangement='FCC')
        pn.add_boundaries()
        B1 = pn.pores("top_boundary")
        B2 = pn.pores("bottom_boundary")
        assert vo.vertex_dimension(pn, B1, B2, 'minmax') == [0.0, 3.0, 0.0,
                                                             3.0, 0.0, 3.0]
        assert pn.num_pores() == 141
