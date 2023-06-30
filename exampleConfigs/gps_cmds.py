def get_monoenergy_5deg20(x,y,z,energy,particle="neutron"):
    """
    Get particle from 5 to 20 degrees. (forward)
    common particles: neutron, pi-, e-, mu-
    common positions:
      z = 870. (back Hcal)
    """
    return [
        "/gps/source/add 1.",
        "/gps/particle %s"%particle,
        "/gps/pos/type Plane",
        "/gps/direction 0 0 1",
        "/gps/pos/shape Square",
        "/gps/pos/centre %.2f %.2f %.2f mm"%(x,y,z),
        "/gps/pos/halfx 1 mm",
        "/gps/pos/halfy 1 mm",
        "/gps/ene/type Mono",
        "/gps/energy %.2f MeV"%energy,
        "/gps/ang/type iso",
        "/gps/ang/mintheta 2.8 rad", # about 20 deg
        "/gps/ang/maxtheta 3.05 rad",
    ]

def get_distenergy_5deg20(x,y,z,energy_min,energy_max,particle="neutron"):
    return [
        "/gps/source/add 1.",
        "/gps/particle %s"%particle,
        "/gps/pos/type Plane",
        "/gps/direction 0 0 1",
        "/gps/pos/shape Square",
        "/gps/pos/centre %.2f %.2f %.2f mm"%(x,y,z),
        "/gps/pos/halfx 1 mm",
        "/gps/pos/halfy 1 mm",
        "/gps/ene/type Lin",
        "/gps/ene/min %.2f GeV"%energy_min,
        "/gps/ene/max %.2f GeV"%energy_max,
        "/gps/ene/gradient 0",
        "/gps/ene/intercept 1",
        "/gps/ang/type iso",
        "/gps/ang/mintheta 2.8 rad", # about 20 deg
        "/gps/ang/maxtheta 3.05 rad",
    ]
