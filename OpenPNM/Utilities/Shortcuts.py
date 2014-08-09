import string, random
import OpenPNM

def solve_linear(pn, ics):
    # circumvent bug with naming by creating random names
    name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
    pseudo = OpenPNM.Fluids.GenericFluid(network=pn, name=name, loglevel=0)
    pseudo.set_throat_data(prop='electronic_conductance', data=1)

    alg = OpenPNM.Algorithms.OhmicConduction(network=pn, name='alg', loglevel=0)
    alg.set_pore_data(prop='BCval', data=ics)
    alg.set_pore_info(label='Dirichlet', locations=ics.nonzero())
    alg.run(active_fluid=pseudo)
    alg.update()

    out = pseudo.get_pore_data(prop='voltage')
    return out