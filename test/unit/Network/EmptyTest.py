import OpenPNM


class EmptyTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Empty(Np=100, Nt=200)

    def test_queries(self):
        assert self.net.Np == 100
        assert self.net.Nt == 200
