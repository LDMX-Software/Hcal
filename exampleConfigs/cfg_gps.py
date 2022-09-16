"""
Pre-defined configurations to test clustering in the back Hcal
"""

from LDMX.Framework import ldmxcfg

p=ldmxcfg.Process("backhcal")
p.libraries.append("libSimCore.so")
p.libraries.append("libHcal.so")
p.libraries.append("libEcal.so")

possible_configs = [
    "gun",
    "mono_gps",
    "gps",
]

import argparse, sys
parser = argparse.ArgumentParser(f'ldmx fire {sys.argv[0]} {sys.argv[1]}')
parser.add_argument('nevents',type=int) # number of events
parser.add_argument('energy',type=float) # energy (MeV)
parser.add_argument('--particle',type=str,
                    choices=['neutron','mu-','pi-','e-'],
                    required=True)
parser.add_argument('--config',type=str,
                    choices=possible_configs,
                    required=True)
parser.add_argument('--x',type=float,
                    default=100.)
parser.add_argument('--y',type=float,
                    default=100.)
parser.add_argument('--energy-max',type=float,
                    default=None)
parser.add_argument('--npart',type=int,
                    default=1)
arg = parser.parse_args()

config = arg.config
energy = arg.energy
particle = arg.particle
energy_max = arg.energy_max
npart = arg.npart

if npart>1 and config=="gun":
    print('No gun allowed with more than one particle')
    exit()
    
from LDMX.SimCore import simulator
from LDMX.SimCore import generators

sim = simulator.simulator("")
sim.setDetector( 'ldmx-det-v12' , True )
sim.runNumber = 0
sim.description = f"{config}"
sim.beamSpotSmear = [20., 80., 0.] #mm
# sim.verbosity = 1

z_back_hcal = 870.
#z_back_hcal = 690.6

from gps_cmds import *
if "gps" in config:
    gpsCmds = []
    if config=="mono_gps":
        gpsCmds = get_monoenergy_5deg20(arg.x,arg.y,z_back_hcal,energy,particle)
        if npart==2:
            gpsCmds += ["/gps/source/multiplevertex true"]
            gpsCmds += get_monoenergy_5deg20(-(arg.x),-(arg.y),z_back_hcal,energy,particle)
    elif config=="gps" and energy_max is not None:
        gpsCmds = get_distenergy_5deg20(arg.x,arg.y,z_back_hcal,energy,energy_max,particle)
        if npart==2:
             gpsCmds += ["/gps/source/multiplevertex true"]
             gpsCmds = get_distenergy_5deg20(-arg.x,-arg.y,z_back_hcal,energy,energy_max,particle)
             
    myGPS = generators.gps( 'myGPS' , gpsCmds )
    print(myGPS)
    sim.generators.append(myGPS)
else:
    particle_gun = generators.gun( "gun_hcal")
    particle_gun.particle = particle
    particle_gun.position = [ 0., 0., z_back_hcal] # mm
    particle_gun.direction = [ 0., 0., 1]
    particle_gun.energy = float(energy/1000.) # GeV
    myGen = particle_gun
    print(myGen)
    sim.generators.append(myGen)

p.outputFiles=['data/%s_%i%s_%.2fmev.root'%(config,npart,particle,energy)]
p.maxEvents = arg.nevents
p.logFrequency = 100
p.sequence=[sim]

import LDMX.Ecal.digi as ecal_digi
import LDMX.Hcal.digi as hcal_digi
from LDMX.Ecal import EcalGeometry
geom = EcalGeometry.EcalGeometryProvider.getInstance()
from LDMX.Hcal import HcalGeometry
geom = HcalGeometry.HcalGeometryProvider.getInstance()
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.ecal_hardcoded_conditions

digit = hcal_digi.HcalDigiProducer()
recon = hcal_digi.HcalRecProducer()

p.sequence.extend([
    ecal_digi.EcalDigiProducer(),
    ecal_digi.EcalRecProducer(),
    digit,
    recon,
])

p.keep = [
    "keep HcalScoringPlaneHits.*",
    "keep EcalScoringPlaneHits.*",
    "drop HcalSimHits.*",
    "drop EcalSimHits.*",
    "drop HcalDigis.*",
    "drop EcalDigis.*",
]
