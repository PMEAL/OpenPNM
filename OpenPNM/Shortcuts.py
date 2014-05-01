import OpenPNM

def solve_linear(pn, ics):
    pseudo = OpenPNM.Fluids.GenericFluid(network=pn, name='none', loglevel=0)
    pseudo.set_throat_data(prop='electronic_conductance', data=1)

    alg = OpenPNM.Algorithms.OhmicConduction(network=pn, name='alg', loglevel=0)
    alg.set_pore_data(prop='BCval', data=ics)
    alg.set_pore_info(label='Dirichlet', locations=ics.nonzero())
    alg.run(active_fluid=pseudo)
    alg.update()

    return pseudo.get_pore_data(prop='voltage')

if __name__ == '__main__':
    pn = OpenPNM.Network.Cubic(name='net')
    pn.generate(add_boundaries=False)

    x,y,z = pn.get_pore_data(prop='coords').T
    left = x==x.min()
    right = x==x.max()
    ics = 2*left + 1*right

    sol = solve_linear(pn, ics)
    print sol.shape