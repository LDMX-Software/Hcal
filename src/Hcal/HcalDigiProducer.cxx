/**
 * @file HcalDigiProducer.cxx
 * @brief Class that performs basic HCal digitization
 * @author Cameron Bravo, SLAC National Accelerator Laboratory
 * @author Tom Eichlersmith, University of Minnesota
 * @author Cristina Suarez, Fermi National Accelerator Laboratory
 * @author Ralf Ehrlich, University of Virginia
 */

#include "Hcal/HcalDigiProducer.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include "Framework/RandomNumberSeedService.h"

namespace hcal {

HcalDigiProducer::HcalDigiProducer(const std::string& name,
                                   framework::Process& process)
    : Producer(name, process),
      engine_(0),  //FIXME: needs a seed
      randFlat_(engine_), randGaussQ_(engine_), randPoissonQ_(engine_) {
  /*
   * Noise generator by default uses a Gausian model for noise
   * i.e. It assumes the noise is distributed around a mean (setPedestal)
   * with a certain RMS (setNoise) and then calculates
   * how many hits should be generated for a given number of empty
   * channels and a minimum readout value (setNoiseThreshold)
   */
  noiseGenerator_ = std::make_unique<ldmx::NoiseGenerator>();
}

void HcalDigiProducer::configure(framework::config::Parameters& ps) {
  // settings of readout chip
  //  used  in actual digitization
  auto hgcrocParams = ps.getParameter<framework::config::Parameters>("hgcroc");
  hgcroc_ = std::make_unique<ldmx::HgcrocEmulator>(hgcrocParams);
  clockCycle_ = hgcrocParams.getParameter<double>("clockCycle");
  nADCs_ = hgcrocParams.getParameter<int>("nADCs");
  iSOI_ = hgcrocParams.getParameter<int>("iSOI");
  noise_ = hgcrocParams.getParameter<bool>("noise");
  photongen_ = hgcrocParams.getParameter<bool>("photongen");

  // collection names
  inputCollName_ = ps.getParameter<std::string>("inputCollName");
  inputPassName_ = ps.getParameter<std::string>("inputPassName");
  digiCollName_ = ps.getParameter<std::string>("digiCollName");

  // physical constants
  //  used to calculate unit conversions
  MeV_ = ps.getParameter<double>("MeV");
  attlength_ = ps.getParameter<double>("attenuationLength");

  // Time -> clock counts conversion
  //  time [ns] * ( 2^10 / max time in ns ) = clock counts
  ns_ = 1024. / clockCycle_;

  // configure scintillation yield
  scintillationYield_= ps.getParameter<double>("scintillationYield");
  scintillationYieldSigma_= ps.getParameter<double>("scintillationYieldSigma");
  scintillationYieldCutoffLow_= ps.getParameter<double>("scintillationYieldCutoffLow");
  scintillationYieldCutoffHigh_= ps.getParameter<double>("scintillationYieldCutoffHigh");

  if(photongen_) {
    // configure photon generator (load lookup table for each scintillator length)
    // reflector type:
    //  0: the counter has SiPMs on both ends, i.e. no reflector or absorber on either side.
    //  1: the counter has a SiPM on one side, and a reflector on the other side.
    //  2: the counter has a SiPM on one side, and an absorber on the other side.
    //  a negative sign indicates that the absorber is at the negative end of the counter.
    std::vector<std::string> lookupTableNames  = ps.getParameter<std::vector<std::string>>("lookupTables");
    for(size_t i=0; i<lookupTableNames.size(); ++i)
    {
      std::shared_ptr<HcalPhotonGenerator> photonGenerator = std::shared_ptr<HcalPhotonGenerator>(new HcalPhotonGenerator(randFlat_, randGaussQ_, randPoissonQ_));
      photonGenerator->LoadLookupTable(lookupTableNames.at(i),1);
      photonGenerator->SetScintillationYield(scintillationYield_);
      int length=(int)(photonGenerator->GetScintillatorLength()+0.5); //rounded to closest int to avoid mismatches due to precission of floating point numbers
      int reflectorType=photonGenerator->GetReflectorType();
      photonGenerators_[std::pair<int,int>(length,reflectorType)] = photonGenerator;
    }
  
    // configure charge generator
    singlePixelPeakVoltage_ = ps.getParameter<double>("singlePixelPeakVoltage"); 
    deadSiPMProbability_ = ps.getParameter<double>("deadSiPMProbability");         //0.01 ???

    // from Mu2e CRV SiPMs (Hamamatsu S13360-2050VE)
    // FIXME: Replace with SiPMS used for Hcal
    int nPixelsX = ps.getParameter<int>("nPixelsX");                               //40
    int nPixelsY = ps.getParameter<int>("nPixelsY");                               //40
    double overvoltage = ps.getParameter<double>("overvoltage");                   //3.0V
    double timeConstant = ps.getParameter<double>("timeConstant");                 //12.0ns
    double capacitance = ps.getParameter<double>("capacitance");                   //8.84e-14F (per pixel)
    double digitizationStart_ = ps.getParameter<double>("digitizationStart");      //0ns
    double digitizationEnd_ = ps.getParameter<double>("digitizationEnd");          //100ns  ???
    //FIXME: How to use vector<pair<int,int>> as parameter?
    std::vector<std::vector<int>> inactivePixelsTmp = ps.getParameter<std::vector<std::vector<int>>>("inactivePixels");  //{18,18},....,{21,21} 
    std::vector<std::pair<int,int>> inactivePixels;
    for(size_t i=0; i<inactivePixelsTmp.size(); ++i) 
      inactivePixels.emplace_back(inactivePixelsTmp[i].at(0),inactivePixelsTmp[i].at(1));
  
    HcalChargeGenerator::ProbabilitiesStruct probabilities;
    probabilities._avalancheProbParam1 = ps.getParameter<double>("AvalancheProbParam1");  //0.65
    probabilities._avalancheProbParam2 = ps.getParameter<double>("AvalancheProbParam2");  //2.7
    probabilities._trapType0Prob = ps.getParameter<double>("TrapType0Prob");              //0
    probabilities._trapType1Prob = ps.getParameter<double>("TrapType1Prob");              //0
    probabilities._trapType0Lifetime = ps.getParameter<double>("TrapType0Lifetime");      //5.0ns
    probabilities._trapType1Lifetime = ps.getParameter<double>("TrapType1Lifetime");      //50.0ns
    probabilities._thermalRate = ps.getParameter<double>("ThermalRate");                  //3.0e-4 ns^-1   300MHz for entire SiPM
    probabilities._crossTalkProb = ps.getParameter<double>("CrossTalkProb");              //0.05

    std::string photonMapFileName = ps.getParameter<std::string>("photonMap");
    chargeGenerator_ = std::shared_ptr<HcalChargeGenerator>(new HcalChargeGenerator(randFlat_, randPoissonQ_, photonMapFileName));
    chargeGenerator_->SetSiPMConstants(nPixelsX, nPixelsY, overvoltage, timeConstant, capacitance, probabilities, inactivePixels);
  }

  // configure generator that will produce noise hits in empty channels
  double readoutThreshold = ps.getParameter<double>("avgReadoutThreshold");
  double gain = ps.getParameter<double>("avgGain");
  double pedestal = ps.getParameter<double>("avgPedestal");
  // rms noise in mV
  noiseGenerator_->setNoise(gain*ps.getParameter<double>("avgNoiseRMS"));
  // mean noise amplitude (if using Gaussian Model for the noise) in mV
  noiseGenerator_->setPedestal(gain * pedestal);
  // threshold for readout in mV
  noiseGenerator_->setNoiseThreshold(gain * readoutThreshold);
}

void HcalDigiProducer::produce(framework::Event& event) {
  // Handle seeding on the first event
  if (!noiseGenerator_->hasSeed()) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    noiseGenerator_->seedGenerator(
        rseed.getSeed("HcalDigiProducer::NoiseGenerator"));
  }
  if (noiseInjector_.get() == nullptr) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    noiseInjector_ = std::make_unique<TRandom3>(
        rseed.getSeed("HcalDigiProducer::NoiseInjector"));
  }
  if (!hgcroc_->hasSeed()) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    hgcroc_->seedGenerator(rseed.getSeed("HcalDigiProducer::HgcrocEmulator"));
  }

  // Get the Hgcroc Conditions
  hgcroc_->condition(
      getCondition<conditions::DoubleTableCondition>("HcalHgcrocConditions"));

  // Get the Hcal Geometry
  const auto& hcalGeometry = getCondition<ldmx::HcalGeometry>(
      ldmx::HcalGeometry::CONDITIONS_OBJECT_NAME);

  // Get the Hcal SimHits
  auto hcalSimHits{event.getCollection<ldmx::SimCalorimeterHit>(inputCollName_, inputPassName_)};
  
  // Get the map of simulated particles
  auto particleByTrackID{event.getMap<int, ldmx::SimParticle>("SimParticles")};

  // Scintillation yield for each HcalID
  std::map<unsigned int, double> scintillationYieldsByID;

  // Photon arrival times for each HcalDigiID 
  std::map<unsigned int, std::vector<double>> photonTimesByID;

  for (int section = 0; section < hcalGeometry.getNumSections(); ++section) {
    for (int layer = 1; layer <= hcalGeometry.getNumLayers(section); ++layer) {
      for(int strip = 0; strip < hcalGeometry.getNumStrips(section, layer); ++strip) {
        auto detID = ldmx::HcalID(section, layer, strip);
	unsigned int hitID = detID.raw();

	// randomly set a scintillation yield around the mean scintillation yield for each scintillator
	double adjustedYield = 0;
	do
	  {
	    adjustedYield = randGaussQ_.fire(scintillationYield_, scintillationYield_*scintillationYieldSigma_);
	  } while(adjustedYield < scintillationYield_*scintillationYieldCutoffLow_ ||
		  adjustedYield > scintillationYield_*scintillationYieldCutoffHigh_);
	scintillationYieldsByID[hitID] = adjustedYield;

	for(int end = 0; end <2; ++end) {
	  photonTimesByID.emplace(4, std::vector<int>());
	}
	  
      }
    }
  }

  // SimHit map for existent HcalIDs
  std::map<unsigned int, std::vector<const ldmx::SimCalorimeterHit*>> hitsByID;
  for (auto const& simHit : hcalSimHits) {
    unsigned int hitID = simHit.getID();
    auto idh = hitsByID.find(hitID);
    if (idh == hitsByID.end()) {
      hitsByID[hitID] = std::vector<const ldmx::SimCalorimeterHit*>(1, &simHit);
    } else {
      idh->second.push_back(&simHit);
    }
  }

  if(photongen_) {
    for (auto const& simBar : hitsByID) {
      unsigned int hitID = simBar.first;
      ldmx::HcalID detID(hitID);
      int section = detID.section();
      int layer = detID.layer();
      int strip = detID.strip();
      ldmx::HcalDigiID posEnd(section, layer, strip, 0);
      ldmx::HcalDigiID negEnd(section, layer, strip, 1);
      unsigned int posEndID = posEnd.raw();
      unsigned int negEndID = negEnd.raw();

      // start a photonGenerator based on scintillator length and reflector type
      int length = (int)(hcalGeometry.getScintillatorLength(detID));
      int reflectorType = (section==ldmx::HcalID::HcalSection::BACK ? 0:2); 
      auto photonGenerator = photonGenerators_.find(std::pair<int,int>(length, reflectorType));
      if(photonGenerator == photonGenerators_.end())
        throw std::runtime_error("HcalDigiProducer::produce: Found an a scintillator for which we don't have a lookup table");
      photonGenerator->second->SetScintillationYield(scintillationYieldsByID[hitID]);

      // loop over all simhits 
      for (auto psimHit : simBar.second) {
	const ldmx::SimCalorimeterHit& simHit = *psimHit;

	int nContributions = simHit.getNumberOfContribs();
	std::cout << " Number of contributions " << nContributions << std::endl;
	///if(nContributions==0) continue;  //There should only be 1 contributions for Hcal hits

        double charge = particleByTrackID[simHit.getContrib(0).trackID].getCharge();

	// get simhit's pre- and post-step positions to the hits in the coordinate frame of the sensitive volume
	// XYZ position: scintillator thickness, width and length
	CLHEP::Hep3Vector pos1Local, pos2Local;
	for(int i=0; i<3; ++i) {
          pos1Local[i] = simHit.getPreStepPosition().at(i);
          pos2Local[i] = simHit.getPostStepPosition().at(i);
        }
	// FIXME: where do 10,24 and length/2 are coming from?
	if(std::fabs(pos1Local[0])>10.0) throw std::runtime_error("step point bigger than lookup table thickness");
        if(std::fabs(pos2Local[0])>10.0) throw std::runtime_error("step point bigger than lookup table thickness");
        if(std::fabs(pos1Local[1])>25.0) throw std::runtime_error("step point bigger than lookup table width");
        if(std::fabs(pos2Local[1])>25.0) throw std::runtime_error("step point bigger than lookup table width");
        if(std::fabs(pos1Local[2])>length/2.0) throw std::runtime_error("step point bigger than lookup table length");
        if(std::fabs(pos2Local[2])>length/2.0) throw std::runtime_error("step point bigger than lookup table length");
		
	photonGenerator->second->MakePhotons(
                                             pos1Local,
                                             pos2Local,
                                             simHit.getPreStepTime(),
                                             simHit.getPostStepTime(),
                                             simHit.getVelocity()/CLHEP::c_light,
                                             charge,
                                             simHit.getEdep(),
                                             simHit.getPathLength(),
                                             reflectorType
                                             );

	// get arrival times of the photon
        // FIXME: What do positive and negative ends mean in the photonGenerator?
        // SiPM1 is at the positive end in the photonGenerator
        const std::vector<double> &timesPosEnd = photonGenerator->second->GetArrivalTimes(1);
        // SiPM0 is at the negative end in the photonGenerator
        const std::vector<double> &timesNegEnd = photonGenerator->second->GetArrivalTimes(0);

	// TODO: implement function in HcalGeometry to deduce where SiPMS are
	switch (section)
          {
          case ldmx::HcalID::HcalSection::BACK   : photonTimesByID[posEndID].insert(photonTimesByID[posEndID].end(),timesPosEnd.begin(),timesPosEnd.end());
            photonTimesByID[negEndID].insert(photonTimesByID[negEndID].end(),timesNegEnd.begin(),timesNegEnd.end());
            break;
          case ldmx::HcalID::HcalSection::TOP    :
          case ldmx::HcalID::HcalSection::RIGHT  : photonTimesByID[negEndID].insert(photonTimesByID[negEndID].end(),timesNegEnd.begin(),timesNegEnd.end());
            break;
          case ldmx::HcalID::HcalSection::BOTTOM :
          case ldmx::HcalID::HcalSection::LEFT   : photonTimesByID[posEndID].insert(photonTimesByID[posEndID].end(),timesPosEnd.begin(),timesPosEnd.end());
            break;
          }
      } // end loop over simHits
    } // end loop over IDs
  } 

  // Empty collection to be filled
  ldmx::HgcrocDigiCollection hcalDigis;
  hcalDigis.setNumSamplesPerDigi(nADCs_);
  hcalDigis.setSampleOfInterestIndex(iSOI_);
  typedef std::vector<ldmx::HgcrocDigiCollection::Sample> digisType;
  std::vector<std::pair<ldmx::HcalDigiID,digisType> > digisToAdd;

  // HGCROC Emulation 
  // loop over all possible HcalDigiIDs in order to also catch SiPMs that don't have hits, but may get dark counts
  for (auto const& simBar : scintillationYieldsByID) {
    ldmx::HcalID detID(simBar.first);
    int section = detID.section();
    int layer = detID.layer();
    int strip = detID.strip();

    // find this HcalDigiID of this SiPM in the photonTimes map
    // and fill the new photon times/index vector
    //auto times = photonTimes.find(channelID);
    //if (times!=photonTimes.end())
  }

  // scintillator orientation (0: along x-axis, 1: along y-axis, 2: along z-axis)
  //const auto orientation{hcalGeometry.getScintillatorOrientation(detID)};
    
  event.add(digiCollName_, hcalDigis);
  
  return;
}  // produce
  
}  // namespace hcal

DECLARE_PRODUCER_NS(hcal, HcalDigiProducer);
