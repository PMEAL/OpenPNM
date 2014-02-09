
from .__GenericFluid__ import GenericFluid

class Air(GenericFluid):
    r"""
    Creates Fluid object with a default name 'air' and preset values
    """
    def __init__(self,**kwargs):
        super(Air,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.set_pore_data(prop='Tc',data=132.65)
        self.set_pore_data(prop='Pc',data=3.771e6)
        self.set_pore_data(prop='MW',data=0.0291)
        self.add_method(prop='diffusivity',model='Fuller',MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
        self.add_method(prop='viscosity',model='Reynolds',uo=0.001,b=0.1)
        self.add_method(prop='molar_density',model='ideal_gas',R=8.314)
        self.regenerate()

if __name__ =="__main__":
    print('no tests yet')