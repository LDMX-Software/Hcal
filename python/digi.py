"""Package to configure the HCal digitization pipeline

All classes are derived versions of LDMX.Framework.ldmxcfg.Producer
with helpful member functions.

Two module-wide parameters are defined.

Attributes
----------
nPEPerMIP: float
    Number of photo-electrons (PEs) created for each MIP 
mipEnergy: float
    Energy [MeV] of a single MIP 
"""

from LDMX.Framework.ldmxcfg import Producer

from LDMX.Tools.HgcrocEmulator import HgcrocEmulator

nPEPerMIP = 68. #PEs created per MIP 
mipEnergy = 4.66 #MeV - measured 1.4 MeV for a 6mm thick tile, so for 20mm bar = 1.4*20/6      

class HcalHgcrocEmulator(HgcrocEmulator) :
    """
    Get an HGCROC emulator and configure for the HCal specifically
    This sets the pulse shape parameters to the ones from a fit
    to a test readout of an HCal module and then thresholds to the
    default construction.
    Noise RMS is calculated using the voltage of 0.02 PEs.
    """

    def __init__(self) :
        super().__init__()

        # SOI
        # Sample of interest (will have double of samples (6) after pulse peak)
        self.iSOI = 3

        # nADCs
        self.nADCs = 10

        # set pulse shape parameters
        self.rateUpSlope = -0.1141
        self.timeUpSlope = -9.897
        self.rateDnSlope = 0.0279
        self.timeDnSlope = 45.037
        self.timePeak    = 12.698 # the time such that with [parameter 4]=0, the pulse peaks at t=0

        # noise (0.02PE)
        self.noiseRMS = self.calculateVoltageHcal(0.02) # mV

    def calculateVoltageHcal(self, PE) :
        """Calculate the voltage signal [mV] of the input number of photo-electrons (PEs)
        Assuming that 1 PE ~ 5mV
        This translates to (68/4.66)*5 = 73 PE/MeV
        Parameters
        ----------
        PE : int
             Number of photo electrons
        """
        return PE*(5/1)
    
class HcalDigiProducer(Producer) :
    """Configuration for HcalDigiProducer

    Attributes
    ----------
    hgcroc : HgcrocEmulator
        Configuration for the chip emulator
    MeV : float
        Conversion between energy [MeV] and voltage [mV]
    inputCollName : str
        Name of input collection  
    inputPassName : str
        Name of input pass 
    digiCollName : str    
        Name of digi collection                                                                                                                                                                          
    """

    def __init__(self, instance_name = 'hcalDigis') :
        super().__init__(instance_name , 'hcal::HcalDigiProducer','Hcal')

        self.hgcroc = HcalHgcrocEmulator()

        #Energy -> Volts converstion
        # energy [MeV] ( 1 MIP / energy per MIP [MeV] ) ( voltage per MIP [mV] / 1 MIP ) = voltage [mV]
        # assuming 1 PEs ~ 5mV ->  self.MeV = 72.961 mV/MeV
        self.MeV = (1./mipEnergy)*self.hgcroc.calculateVoltageHcal( nPEPerMIP )

        # attenuation length
        self.attenuationLength = 5.; # in m   

        # scintillation yield
        self.scintillationYield = 30000.; # in photons/MeV

        # lookup tables
        self.lookupTables = [ "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_2000_0",
                              "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_1800_2",
                              "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_1600_2",
                              "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_1200_2",
                              "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_1400_2", 
                              "@CMAKE_INSTALL_PREFIX@/data/Hcal/LookupTable_600_2" ];

        # SiPM single-pixel response
        self.singlePixelPeakVoltage = 5.0; # 5 mV

        # SiPM photon map
        self.photonMap = "@CMAKE_INSTALL_PREFIX@/data/Hcal/photonMap18.root";

        # other SiPM constants
        self.digitizationStart            =  400.0;     # 0ns  ????
        self.digitizationEnd              = 1750.0;     # 100ns  ????
        self.deadSiPMProbability          = 0.01;       # ????
        self.nPixelsX                     = 40;
        self.nPixelsY                     = 40;
        self.inactivePixels               = [ [18,18], [18,19], [18,20], [18,21],
                                              [19,18], [19,19], [19,20], [19,21],
                                              [20,18], [20,19], [20,20], [20,21],
                                              [21,18], [21,19], [21,20], [21,21] ];
        self.overvoltage                  = 3.0;        #V
        self.timeConstant                 = 13.3;       #ns  according to an Hamamatsu example with R_q=150kOhm --> tau=R_q*C=13.3ns
        self.capacitance                  = 8.84e-14;   #F   capacitance of one pixel according to specs

        self.AvalancheProbParam1          = 0.607;      # = p1
        self.AvalancheProbParam2          = 2.7;        # = p2
                                                        # Avalanche probability at over voltage v: p1*(1 - exp(-v/p2))

        self.TrapType0Prob                = 0.0;        # 0.14 (Paul's number)  ????
        self.TrapType1Prob                = 0.0;        # 0.06 (Paul's number)  ????
        self.TrapType0Lifetime            = 5.0;        # ns  ????
        self.TrapType1Lifetime            = 50.0;       # ns  ????

        self.ThermalRate                  = 3.0e-4;     # ns^-1     0.3MHz for entire SiPM
        self.CrossTalkProb                = 0.05;       #

        # avg parameters
        self.avgReadoutThreshold = 4. #ADCs - noise config only
        self.avgGain = 1.2 #noise config only
        self.avgPedestal = 1. #noise config only
        
        # input and output collection name parameters
        self.inputCollName = 'HcalSimHits'
        self.inputPassName = ''
        self.digiCollName = 'HcalDigis'

class HcalRecProducer(Producer) :
    """Configuration for the HcalRecProducer

    Attributes
    ----------
    voltage_per_mip: float
        Conversion from voltage [mV] to number of MIPs
    mip_energy : float
        Copied from module-wide mipEnergy [MeV]
    clock_cycle : float
        Time for one DAQ clock cycle to pass [ns]
    digiCollName : str
        Name of digi collection
    digiPassName : str
        Name of digi pass
    simHitCollName : str
        Name of simHit collection
    simHitPassName : str 
        Name of simHit pass 
    recHitCollName : str
        Name of recHit collection
    """

    def __init__(self, instance_name = 'hcalRecon') : 
        super().__init__(instance_name , 'hcal::HcalRecProducer','Hcal')

        hgcroc = HcalHgcrocEmulator()

        self.voltage_per_mip = (5/1)*(nPEPerMIP) # 5*68 mV/ MIP
        self.mip_energy = mipEnergy #MeV / MIP
        self.clock_cycle = 25. #ns - needs to match the setting on the chip   
        self.pe_per_mip = nPEPerMIP
        
	# attenuation length
        self.attenuationLength = 5.; # in m  
        
        self.digiCollName = 'HcalDigis'
        self.digiPassName = ''
        self.simHitCollName = 'HcalSimHits'
        self.simHitPassName = ''
        self.recHitCollName = 'HcalRecHits'

        # hgcroc parameters:
        self.rateUpSlope = hgcroc.rateUpSlope
        self.timeUpSlope = hgcroc.timeUpSlope
        self.rateDnSlope = hgcroc.rateDnSlope
        self.timeDnSlope = hgcroc.timeDnSlope
        self.timePeak    = hgcroc.timePeak
        self.nADCs       = hgcroc.nADCs

        # avg parameters
        self.avgToaThreshold = 1.6 # mV - correction config only
        self.avgGain = 1.2 # correction config only 
        self.avgPedestal = 1. #noise config only   

class HcalSingleEndRecProducer(Producer) :
    """ Configuration for the single ended Hcal Rec Producer

    Attributes
    ----------
    -  mip_energy : float
       Copied from module-wide mipEnergy [MeV]
    -  clock_cycle : float
       Time for one DAQ clock cycle to pass [ns]
    -  pe_per_mip: float
       number of photo-electrons per MIP
    -  pass_name: str
       Name of digi pass
    -  coll_name: str
       Name of digi collection
    -  rec_coll_name: str
       Name of rechit collection
    """

    def __init__(self, instance_name = 'hcalRecon', pass_name = '', coll_name = 'HcalDigis', rec_coll_name = 'HcalRecHits', rec_pass_name = '') :
        super().__init__(instance_name , 'hcal::HcalSingleEndRecProducer','Hcal')

        self.mip_energy = mipEnergy
        self.clock_cycle = 25.
        self.pe_per_mip = nPEPerMIP
        
        self.coll_name = coll_name
        self.pass_name = pass_name
        self.rec_coll_name = rec_coll_name
        self.rec_pass_name = rec_pass_name

class HcalDoubleEndRecProducer(Producer) :
    """ Configuration for the double ended Hcal Rec Producer
    
    Attributes
    ----------
    -  mip_energy : float
       Copied from module-wide mipEnergy [MeV]
    -  clock_cycle : float
       Time for one DAQ clock cycle to pass [ns]
    -  pe_per_mip: float
       number of photo-electrons per MIP
    -  pass_name: str
       Name of digi pass
    -  coll_name: str
       Name of digi collection
    -  rec_coll_name: str
       Name of rechit collection
    """

    def __init__(self, instance_name = 'hcalDoubleRecon', pass_name = '', coll_name = 'HcalRecHits', rec_coll_name = 'HcalDoubleEndRecHits', rec_pass_name = '') :
        super().__init__(instance_name , 'hcal::HcalDoubleEndRecProducer','Hcal')

        self.mip_energy = mipEnergy
        self.clock_cycle = 25.
        self.pe_per_mip = nPEPerMIP

        self.coll_name = coll_name
        self.pass_name = pass_name
        self.rec_coll_name = rec_coll_name
        self.rec_pass_name = rec_pass_name
