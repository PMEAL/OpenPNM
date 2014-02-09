
from .__GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    Creates Fluid object with a default name 'water' and preset values
    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.set_pore_data(prop='Tc',data=647.096)
        self.set_pore_data(prop='Pc',data=22.06e6)
        self.set_pore_data(prop='MW',data=0.0291)
        self.add_method(prop='diffusivity',model='constant',value=1e-12)
        self.add_method(prop='viscosity',model='constant',value=0.001)
        self.add_method(prop='molar_density',model='constant',value=44445)
        self.add_method(prop='surface_tension',model='constant',value=0.072)
        self.add_method(prop='contact_angle',model='constant',value=110)
        self.regenerate()

if __name__ =="__main__":
    print('no tests yet')
