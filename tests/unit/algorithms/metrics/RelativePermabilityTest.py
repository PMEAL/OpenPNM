import openpnm as op
import scipy as sp
mgr = op.Workspace()


class RelativePermeabilityTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10], spacing=0.0006)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.non_wet_phase = op.phases.Air(network=self.net)
        self.wet_phase = op.phases.Water(network=self.net)
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.non_wet_phase.add_model(propname='throat.hydraulic_conductance',
              model=mod)
        self.wet_phase.add_model(propname='throat.hydraulic_conductance',
              model=mod)
        self.inlet_pores=self.net.pores('left')
        ip = op.algorithms.InvasionPercolation(network=self.net, phase=self.non_wet_phase)
        ip.set_inlets(pores=self.inlet_pores)
        ip.run()
        self.non_wet_phase.update(ip.results())
        

    def test_one_phase_definition(self):
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                 invasion_sequence='invasion_sequence')
        #rp.set_inlets(pores=Finlets_init['x'])
        #rp.set_outlets(pores=Foutlets_init['x'])
        rp.run()
        results=rp.get_Kr_data()
        assert results['results']['krw']['x']!=[]


    def test_overwriting_boundary_faces(self):
        inlets={'x':'left','y':'left','z':'left'}
        outlets={'x':'right','y':'right','z':'right'}
        rp = op.algorithms.metrics.RelativePermeability(network=self.net)
        rp.setup(invading_phase=self.non_wet_phase,
                 defending_phase=self.wet_phase,
                 invasion_sequence='invasion_sequence',
                 flow_inlets=inlets,
                 flow_outlets=outlets)
        rp.run()
        results=rp.get_Kr_data()
        assert results['results']['krw']['x']!=results['results']['krw']['y']
        
        
    def test_lacking_boundary_faces(self):
    inlets={'x':'left'}
    outlets={'x':'right'}
    rp = op.algorithms.metrics.RelativePermeability(network=self.net)
    rp.setup(invading_phase=self.non_wet_phase,
             defending_phase=self.wet_phase,
             invasion_sequence='invasion_sequence',
             flow_inlets=inlets,
             flow_outlets=outlets)
    rp.run()
    results=rp.get_Kr_data()
    assert results['results']['krw']['y']==[]
#    def test_given_length(self):
#        FF = op.algorithms.metrics.FormationFactor(network=self.net)
#        FF.run()
#        val_1 = FF.results['y']
#        FF.set_length(direction='y', length=15*0.0005)
#        FF.run()
#        val_2 = FF.results['y']
#        assert val_1 != val_2
#
#    def test_setting_inlets(self):
#        FF = op.algorithms.metrics.FormationFactor(network=self.net)
#        Ps = self.net.pores('left')
#        self.net.set_label(pores=Ps, label='blah')
#        FF.set_inlets(direction='x', label='blah')
#        FF.run()
#        val_1 = FF.results['x']
#        FF.set_inlets(direction='x', label='left')
#        FF.run()
#        val_2 = FF.results['x']
#        assert val_1 == val_2


if __name__ == '__main__':

    t = RelativePermeabilityTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
