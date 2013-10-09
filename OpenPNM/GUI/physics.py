import OpenPNM

def CapillaryPressureWashburn(network, sigma=0.072, theta=120):
  network.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure().Washburn(**locals())
  return {'network': network}