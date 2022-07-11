import openpnm as op
from numpy.testing import assert_allclose
import chemicals


class LatentHeatTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            chemicals.phase_change.MK,
            chemicals.phase_change.SMK,
            chemicals.phase_change.Velasco,
            chemicals.phase_change.Clapeyron,
            chemicals.phase_change.Riedel,
            chemicals.phase_change.Chen,
            chemicals.phase_change.Vetere,
            chemicals.phase_change.Liu,
            # chemicals.phase_change.Watson,  # Needs Hvap_ref
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_pure_prop_wrapper_wrapper(target=h2o, f=f).mean())
        assert_allclose(vals, 41276, rtol=0.5)


if __name__ == '__main__':

    t = LatentHeatTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
