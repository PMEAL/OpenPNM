import string, random
import OpenPNM.Phases

def solve_linear(pn, ics):
    # circumvent bug with naming by creating random names
    name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
    pseudo = OpenPNM.Phases.GenericPhase(network=pn, name=name, loglevel=0)
    pseudo['throat.electrical_conductance']=1

    alg = OpenPNM.Algorithms.OhmicConduction(network=pn)
    alg['pore.BCval']=ics
    alg['pore.Dirichlet']=ics!=0
    alg.run(active_phase=pseudo)
    alg.return_results()

    out = pseudo['pore.voltage']
    return out