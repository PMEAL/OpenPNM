def only_option(network,
  capillary_pressure=['Washburn', 'Purcell', 'Imbibition'],
  heat_conduction=['None', 'get_pore_temp'],
  mass_transfer=['None', 'conductivity', 'conduit_diffusion'],
  multiphase=['Phase calculator', 'Henry\'s Law', 'Raoult\'s Law'],):

  return {'modified_network':network}