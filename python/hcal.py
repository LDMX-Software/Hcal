"""Configuration for Hcal pipeline"""

from LDMX.Framework import ldmxcfg

class HcalVetoProcessor(ldmxcfg.Producer) :
    """Configuration for veto in HCal

    Sets all parameters to reasonable defaults.

    Examples
    --------
        from LDMX.EventProc.hcal import HcalVetoProcessor
        p.sequence.append( HcalVetoProcessor() )
    """

    def __init__(self,name = 'hcalVeto') :
        super().__init__(name,'hcal::HcalVetoProcessor','Hcal')

        self.pe_threshold = 5.0
        self.max_time = 50.0
        self.max_depth = 4000.0
        self.back_min_pe = 1.


class HcalOldDigiProducer(ldmxcfg.Producer) :
    """Configuration for Digitization producer in the HCal
        Sets all parameters to reasonable defaults.
    Examples
    --------
        from LDMX.EventProc.hcal import HcalDigiProducer
        p.sequence.append( HcalDigiProducer() )
    """

    def __init__(self,name = 'hcalOldDigis') :
        super().__init__(name,'hcal::HcalOldDigiProducer','Hcal')

        self.meanNoise = 0.02
        self.readoutThreshold= 1
        self.strips_side_lr_per_layer = 12
        self.num_side_lr_hcal_layers = 26
        self.strips_side_tb_per_layer = 12
        self.num_side_tb_hcal_layers = 28
        self.strips_back_per_layer = 60 # n strips correspond to 5 cm wide bars
        #self.num_back_hcal_layers = 96
        self.num_back_hcal_layers = 150
        self.super_strip_size = 1 # 1 = 5 cm readout, 2 = 10 cm readout, ...
        self.mev_per_mip = 4.66  # measured 1.4 MeV for a 6mm thick tile, so for 20mm bar = 1.4*20/6
        self.pe_per_mip = 68. # PEs per MIP at 1m (assume 80% attentuation of 1m)
        self.strip_attenuation_length = 5. # this is in m
        self.strip_position_resolution = 150. # this is in mm
        self.sim_hit_pass_name = '' #use any pass available

class HcalMIPTracking(ldmxcfg.Producer) :
    """Configuration for MIP Tracking in HCal

    Sets all parameters to reasonable defaults.

    Examples
    --------
        from LDMX.EventProc.hcal import HcalMIPTracking
        p.sequence.append( HcalMIPTracking() )
    """

    def __init__(self,name = 'hcalMIPTracks') :
        super().__init__(name,'hcal::HcalMIPTracking','Hcal')

        self.BACK_HCAL_START_Z_ = 840.
        self.MIP_MIN_PE_ = 36.
        #self.MIP_MAX_PE_ = 900.
        self.MIP_MAX_PE_ = 1800.
        self.MIN_TRACK_HITS_ = 4
        self.MIN_SEED_HITS_ = 3
        self.MAX_SEED_HIT_ERROR_ = 300.
        self.MAX_TRACK_EXTRAP_SIGMA_ = 5.
        self.MAX_LAYERS_CONSEC_MISSED_ = 2
        self.USE_ISOLATED_HITS_ = False
        self.NUM_HITS_REQ_ = 3
        self.NUM_HITS_IN_GROUP_ = 4
        self.NUM_GROUPS_REQ_ = 4
        self.NUM_GROUPS_PER_LAY_ = 5
        self.NUM_LAY_PER_GROUP_ = 2

        digiprod = HcalDigiProducer()

        self.strips_side_lr_per_layer = digiprod.strips_side_lr_per_layer
        self.num_side_lr_hcal_layers = digiprod.num_side_lr_hcal_layers
        self.strips_side_tb_per_layer = digiprod.strips_side_tb_per_layer
        self.num_side_tb_hcal_layers = digiprod.num_side_tb_hcal_layers
        self.strips_back_per_layer = digiprod.strips_back_per_layer
        self.num_back_hcal_layers = digiprod.num_back_hcal_layers
        self.strip_position_resolution = digiprod.strip_position_resolution
