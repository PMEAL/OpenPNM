import openpnm as op
from numpy.testing import assert_approx_equal, assert_array_almost_equal, assert_allclose
import chemicals


class DensityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325.0  # Pa
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.molar_density'] = 55539.0  # mol/m3
        self.phase['pore.salinity'] = 0.0  # ppt

    def test_standard(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.standard)
        assert_approx_equal(self.phase['pore.density'].mean(), 999.702)

    def test_ideal_gas(self):
        # Water vapor
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.ideal_gas)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 0.73610248)

    def test_water(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.water)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 996.9522)

    def test_generic_chemicals_for_pure_gas_molar_volume(self):
        mods = [
            chemicals.virial.BVirial_Pitzer_Curl,
            chemicals.virial.BVirial_Abbott,
            chemicals.virial.BVirial_Tsonopoulos,
            chemicals.virial.BVirial_Tsonopoulos_extended,
        ]
        n2 = op.phase.Species(network=self.net, species='nitrogen')
        n2['pore.temperature'] = 400
        Vm = []
        for f in mods:
            Vm.append(op.models.phase.chemicals_pure_prop(target=n2, f=f).mean())
        assert_allclose(Vm, 8.795e-6, rtol=.3)

    def test_generic_chemicals_for_pure_liq_molar_volume(self):
        mods = [
            chemicals.volume.Yen_Woods_saturation,
            chemicals.volume.Rackett,
            chemicals.volume.Yamada_Gunn,
            chemicals.volume.Townsend_Hales,
            chemicals.volume.Bhirud_normal,
            chemicals.volume.COSTALD,
            chemicals.volume.Campbell_Thodos,
            # chemicals.volume.SNM0,  # numba version not working
            # chemicals.volume.CRC_inorganic,  # requires rho
            # chemicals.volume.COSTALD_compressed,  # requires Psat
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        Vm = []
        for f in mods:
            Vm.append(op.models.phase.chemicals_pure_prop(target=h2o, f=f).mean())
        assert_allclose(Vm, 1.88e-5, rtol=0.2)


if __name__ == '__main__':

    t = DensityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
