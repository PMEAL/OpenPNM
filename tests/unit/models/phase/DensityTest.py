import openpnm as op
from numpy.testing import assert_approx_equal, assert_allclose
from openpnm.utils import get_mixture_model_args
import chemicals
from thermo import Chemical


class DensityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Species(network=self.net, species='h2o')
        self.phase['pore.salinity'] = 0.0  # ppt

    def test_standard(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.liquid_pure_COSTALD)
        assert_approx_equal(self.phase['pore.density'].mean(), 992.345519756)

    def test_ideal_gas(self):
        # Water vapor
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.ideal_gas)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 0.7367280065145)

    def test_water(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.water_correlation)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 996.9522)

    def test_chemicals_for_pure_gas_molar_volume(self):
        mods = [
            # numba version not working for any
            chemicals.virial.BVirial_Pitzer_Curl,
            chemicals.virial.BVirial_Abbott,
            chemicals.virial.BVirial_Tsonopoulos,
            chemicals.virial.BVirial_Tsonopoulos_extended,
        ]
        n2 = op.phase.Species(network=self.net, species='nitrogen')
        n2['pore.temperature'] = 400
        Vm = []
        for f in mods:
            Vm.append(op.models.phase.chemicals_wrapper(n2, f=f).mean())
        assert_allclose(Vm, 8.795e-6, rtol=.3)

    def test_chemicals_wrapper_for_pure_liq_molar_volume(self):
        mods = [
            chemicals.volume.Yen_Woods_saturation,
            chemicals.volume.Rackett,
            chemicals.volume.Yamada_Gunn,
            chemicals.volume.Townsend_Hales,
            chemicals.volume.Bhirud_normal,
            chemicals.volume.COSTALD,
            # chemicals.volume.Campbell_Thodos,  # numba version not working
            # chemicals.volume.SNM0,  # numba version not working
            # chemicals.volume.CRC_inorganic,  # requires rho
            # chemicals.volume.COSTALD_compressed,  # requires Psat
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        Vm = []
        for f in mods:
            Vm.append(op.models.phase.chemicals_wrapper(h2o, f=f).mean())
        assert_allclose(Vm, 1.88e-5, rtol=0.2)

    def test_chemicals_wrapper_for_pure_liq_with_args(self):
        h2o = op.phase.Species(network=self.net, species='water')
        # Using kwargs to map args to custom propnames
        temp = Chemical('h2o')
        h2o['pore.density'] = temp.rhol
        Vm = op.models.phase.chemicals_wrapper(
            phase=h2o,
            f=chemicals.volume.CRC_inorganic,
            rho0='pore.density',
            k=1,
        )
        assert_allclose(Vm, 1.85309071e-05, rtol=1e-4)
        # Put args directly in target
        h2o['pore.Psat'] = temp.Psat
        h2o['pore.Vs'] = temp.Vms
        Vm = op.models.phase.chemicals_wrapper(
            phase=h2o,
            f=chemicals.volume.COSTALD_compressed,
            rho='pore.density',
        )
        assert_allclose(Vm, 1.61982081e-05, rtol=1e-4)

    def test_chemicals_wrapper_for_liquid_mixture(self):
        h2o = op.phase.Species(network=self.net, species='h2o')
        etoh = op.phase.Species(network=self.net, species='ethanol')
        vodka = op.phase.LiquidMixture(network=self.net, components=[h2o, etoh])
        vodka.x(h2o.name, 0.60)
        vodka.x(etoh.name, 0.40)
        Vm = op.models.phase.chemicals_wrapper(
            phase=vodka,
            f=chemicals.volume.COSTALD_mixture,
        )
        h2o['param.Zr'] = 0.001
        etoh['param.Zr'] = 0.001
        Vm = op.models.phase.chemicals_wrapper(
            phase=vodka,
            f=chemicals.volume.Rackett_mixture,
            Zrs='param.Zr',
        )

    def test_liquid_pure_and_mixture(self):
        pn = op.network.Demo()
        h2o = op.phase.Species(network=pn, species='water')
        h2o.add_model(propname='pore.density',
                      model=op.models.phase.density.liquid_pure_COSTALD)
        Vm = chemicals.COSTALD(
            T=h2o['pore.temperature'][0],
            Tc=h2o['param.critical_temperature'],
            Vc=h2o['param.critical_volume'],
            omega=h2o['param.acentric_factor'],
        )
        rho_ref = chemicals.Vm_to_rho(Vm, h2o['param.molecular_weight'])
        rho_calc = h2o['pore.density'][0]
        assert_allclose(rho_ref, rho_calc, rtol=1e-10, atol=0)

        etoh = op.phase.Species(network=pn, species='ethanol')
        etoh.add_model(propname='pore.density',
                       model=op.models.phase.density.liquid_pure_COSTALD)

        vodka = op.phase.LiquidMixture(network=pn, components=[h2o, etoh])
        vodka.x(h2o.name, 0.5)
        vodka.x(etoh.name, 0.5)
        vodka.add_model(propname='pore.density',
                        model=op.models.phase.density.liquid_mixture_COSTALD)
        args = get_mixture_model_args(
            phase=vodka,
            composition='xs',
            args={
                'Tcs': 'param.critical_temperature',
                'Vcs': 'param.critical_volume',
                'omegas': 'param.acentric_factor',
            })
        Vm = chemicals.COSTALD_mixture(T=vodka['pore.temperature'][0], **args)
        rho_ref = chemicals.Vm_to_rho(Vm, vodka.get_mix_vals('param.molecular_weight')[0])
        rho_calc = vodka['pore.density'][0]
        assert_allclose(rho_ref, rho_calc, rtol=1e-10, atol=0)


if __name__ == '__main__':

    t = DensityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
