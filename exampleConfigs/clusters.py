from LDMX.Framework import ldmxcfg
p=ldmxcfg.Process("neutron")
p.libraries.append("libSimCore.so")
p.libraries.append("libHcal.so")
p.libraries.append("libEcal.so")

import argparse, sys
parser = argparse.ArgumentParser(f'ldmx fire {sys.argv[0]}')
parser.add_argument('input_file',type=str)
arg = parser.parse_args()

p.inputFiles=[arg.input_file]
p.outputFiles=[arg.input_file.replace('.root','_cluster.root')]
# p.maxEvents = 1
p.logFrequency = 100

from LDMX.Hcal import HcalGeometry
geom = HcalGeometry.HcalGeometryProvider.getInstance()

import LDMX.Hcal.cluster as hcal_cluster
clusters = hcal_cluster.HcalNewClusterProducer()
clusters.seed_threshold_2d = 0.1 # MeV
clusters.neighbor_threshold_2d = 0.01 # MeV
clusters.seed_threshold_3d  = 4. # MeV
clusters.neighbor_threshold_3d = 0.2 # MeV
clusters.num_neighbors = 4 # number of neighboring strips
clusters.max_xy_2d = 100. # mm 
clusters.max_xy_3d = 1000. # mm
clusters.max_xy_2d_merge = 100. #mm
clusters.use_toa = True
clusters.layer_parity = 1 # if its 1, odd layers

p.sequence = [clusters]
