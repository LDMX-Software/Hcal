"""Pedestal file extraction, currently for testbeam"""

import argparse, sys

parser = argparse.ArgumentParser(f'ldmx fire {sys.argv[0]}')

parser.add_argument('input_file')
parser.add_argument('--input_file_dpm1',default=None,help="input pedestal raw file for second polarfile")
parser.add_argument('--make_histos',action='store_true',help="Make channel-by-channel histogram files")
parser.add_argument('--calib_file',required=True,help="Calibration output file")
parser.add_argument('--max_events',default=10000,type=int)

arg = parser.parse_args()

from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('hcal_ana')
p.maxEvents = arg.max_events
p.termLogLevel = 0
p.logFrequency = 10

import LDMX.Hcal.hgcrocFormat as hcal_format
import LDMX.Hcal.hcal_ana as hcal_ana
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
#from LDMX.DQM import dqm
from LDMX.Packing import rawio

import os
base_name = os.path.basename(arg.input_file).replace('.raw','')

p.outputFiles = [f'/dev/null']

p.histogramFile = f'hist_{base_name}.root'
tbl = f'{os.environ["LDMX_BASE"]}/ldmx-sw/Hcal/data/testbeam_connections.csv'

digi_name_dpm0 = 'HcalRawDigis'
comments = "Source file "+arg.input_file
if arg.input_file_dpm1:
    digi_name_dpm0 = 'polarfire0Digis'
    comments += " " +arg.input_file_dpm1
    
# sequence
#   1. decode event packet into digi collection
#   2. merge event packets if needed
#   3. process pedestals
p.sequence = [
        hcal_format.HcalRawDecoder(
            input_file = arg.input_file,
            connections_table = tbl,
            output_name = digi_name_dpm0,
            detector_name = 'ldmx-hcal-prototype-v1.0'
        )
]
if arg.input_file_dpm1:
    p.sequence.extend([
        hcal_format.HcalRawDecoder(
            input_file = arg.input_file_dpm1,
            connections_table = tbl,
            output_name = 'polarfire1Digis',
            detector_name = 'ldmx-hcal-prototype-v1.0'
        ),
        hcal_format.HcalAlignPolarfires(
            input_names = ['polarfire0Digis','polarfire1Digis'],
            output_name = 'HcalRawDigis'
        ),
        #dqm.NtuplizeHgcrocDigiCollection(
        #    input_name = "HcalRawDigis",
        #    already_aligned = True
        #)
    ])
    
p.sequence.append(
    hcal_ana.HcalPedestalAnalyzer(
        input_name = 'HcalRawDigis',
        output_file = arg.calib_file,
        make_histos = arg.make_histos,
        comments = "Source file "+arg.input_file
    )
)
