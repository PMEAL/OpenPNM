import OpenPNM

def CapillaryPressureWashburn(net, sigma=0.072, theta=120):
  OpenPNM.Physics.CapillaryPressure.Washburn(**locals())
  return {'net': net}