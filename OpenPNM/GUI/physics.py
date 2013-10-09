import OpenPNM

def CapillaryPressureWashburn(net, sigma=0.072, theta=120):
  net.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure().Washburn(**locals())
  return {'net': net}