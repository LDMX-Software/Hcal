/**
 * @file HcalDigiProducer.cxx
 * @brief Class that performs basic HCal digitization
 * @author Cameron Bravo, SLAC National Accelerator Laboratory
 * @author Tom Eichlersmith, University of Minnesota
 * @author Cristina Suarez, Fermi National Accelerator Laboratory
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

  // configure photon generator
  scintillationYield_= ps.getParameter<double>("scintillationYield");
  scintillationYieldSigma_= ps.getParameter<double>("scintillationYieldSigma");
  scintillationYieldCutoffLow_= ps.getParameter<double>("scintillationYieldCutoffLow");
  scintillationYieldCutoffHigh_= ps.getParameter<double>("scintillationYieldCutoffHigh");
  std::vector<std::string> lookupTableNames  = ps.getParameter<std::vector<std::string>>("lookupTables");

  // load lookup table for each scintillator length
  for(size_t i=0; i<lookupTableNames.size(); ++i)
  {
    std::shared_ptr<HcalPhotonGenerator> photonGenerator = std::shared_ptr<HcalPhotonGenerator>(new HcalPhotonGenerator(randFlat_, randGaussQ_, randPoissonQ_));
    photonGenerator->LoadLookupTable(lookupTableNames.at(i),1);
    photonGenerator->SetScintillationYield(scintillationYield_);
    int length=(int)(photonGenerator->GetScintillatorLength()+0.5); //rounded to closest int to avoid mismatches due to precission of floating point numbers
    int reflectorType=photonGenerator->GetReflectorType();
    photonGenerators_[std::pair<int,int>(length,reflectorType)]=photonGenerator;
  }
  scintillationYieldsAdjusted_.clear();

  // Configure charge generator
  singlePixelPeakVoltage_ = ps.getParameter<double>("singlePixelPeakVoltage");
  deadSiPMProbability_ = ps.getParameter<double>("deadSiPMProbability");         //0.01 ???
  int nPixelsX = ps.getParameter<int>("nPixelsX");                               //40
  int nPixelsY = ps.getParameter<int>("nPixelsY");                               //40
  double overvoltage = ps.getParameter<double>("overvoltage");                   //3.0V
  double timeConstant = ps.getParameter<double>("timeConstant");                 //12.0ns
  double capacitance = ps.getParameter<double>("capacitance");                   //8.84e-14F (per pixel)
  digitizationStart_ = ps.getParameter<double>("digitizationStart");             //0ns
  digitizationEnd_ = ps.getParameter<double>("digitizationEnd");                 //100ns  ???
  //FIXME: How to use vector<pair<int,int>> as parameter?
  std::vector<std::vector<int>> inactivePixelsTmp = ps.getParameter<std::vector<std::vector<int>>>("inactivePixels");  //{18,18},....,{21,21} 
  std::vector<std::pair<int,int>> inactivePixels;
  for(size_t i=0; i<inactivePixelsTmp.size(); ++i) inactivePixels.emplace_back(inactivePixelsTmp[i].at(0),inactivePixelsTmp[i].at(1));
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

  // Configure generator that will produce noise hits in empty channels
  double readoutThreshold = ps.getParameter<double>("avgReadoutThreshold");
  double gain = ps.getParameter<double>("avgGain");
  double pedestal = ps.getParameter<double>("avgPedestal");
  // rms noise in mV
  noiseGenerator_->setNoise(
      hgcrocParams.getParameter<double>("noiseRMS"));  // rms noise in mV
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

  // Empty collection to be filled
  ldmx::HgcrocDigiCollection hcalDigis;
  hcalDigis.setNumSamplesPerDigi(nADCs_);
  hcalDigis.setSampleOfInterestIndex(iSOI_);

  std::map<unsigned int, std::vector<const ldmx::SimCalorimeterHit*>> hitsByID;  //old code
  std::map<ldmx::HcalDigiID, std::vector<double>> photonTimes;  //new code

  // get simulated hcal hits from Geant4 and group them by id
  auto hcalSimHits{event.getCollection<ldmx::SimCalorimeterHit>(inputCollName_, inputPassName_)};

  // get map of simulated particles
  auto particleMap{event.getMap<int, ldmx::SimParticle>("SimParticles")};

  for (auto const& simHit : hcalSimHits) {
    // get ID
    unsigned int hitID = simHit.getID();

#define NEWCODE
#ifdef NEWCODE

    //******************//
    //***  new code  ***//
    //******************//

    ldmx::HcalID detID(simHit.getID());
    int section = detID.section();
    int layer = detID.layer();
    int strip = detID.strip();
    ldmx::HcalDigiID posEndID(section, layer, strip, 0);
    ldmx::HcalDigiID negEndID(section, layer, strip, 1);

    //TODO: temporary implementation for getting the scintillator length
    int length=(int)(2.0*hcalGeometry.getHalfTotalWidth(section,(layer%2==1?layer+1:layer-1)) + 0.5); //rounded to closest int to avoid mismatches due to precission of floating point numbers
    int reflectorType=(section==ldmx::HcalID::HcalSection::BACK?0:2);
    auto photonGenerator = photonGenerators_.find(std::pair<int,int>(length, reflectorType));
    if(photonGenerator==photonGenerators_.end()) throw std::runtime_error("HcalDigiProducer::produce: Found an a scintillator for which we don't have a lookup table");

    //randomly set a scintillation yield around the mean scintillation yield for each scintillator
    //needs to be done here, because the geometry is not available in configure()
    auto currentScintillationYield = scintillationYieldsAdjusted_.find(detID);
    if(currentScintillationYield==scintillationYieldsAdjusted_.end())
    {
      double adjustedYield=0;
      do
      {
        adjustedYield=randGaussQ_.fire(scintillationYield_, scintillationYield_*scintillationYieldSigma_);
      } while(adjustedYield<scintillationYield_*scintillationYieldCutoffLow_ ||
              adjustedYield>scintillationYield_*scintillationYieldCutoffHigh_);

      currentScintillationYield = scintillationYieldsAdjusted_.emplace(detID,adjustedYield).first;
    }
    photonGenerator->second->SetScintillationYield(currentScintillationYield->second);

    CLHEP::Hep3Vector pos1Local, pos2Local;
    for(int i=0; i<3; ++i)
    {
      pos1Local[i]=simHit.getPreStepPosition().at(i);
      pos2Local[i]=simHit.getPostStepPosition().at(i);
  if(std::fabs(pos1Local[0])>10.0) throw std::runtime_error("step point bigger than lookup table thickness");
  if(std::fabs(pos2Local[0])>10.0) throw std::runtime_error("step point bigger than lookup table thickness");
  if(std::fabs(pos1Local[1])>25.0) throw std::runtime_error("step point bigger than lookup table width");
  if(std::fabs(pos2Local[1])>25.0) throw std::runtime_error("step point bigger than lookup table width");
  if(std::fabs(pos1Local[2])>length/2.0) throw std::runtime_error("step point bigger than lookup table length");
  if(std::fabs(pos2Local[2])>length/2.0) throw std::runtime_error("step point bigger than lookup table length");
    }

    int nContributions = simHit.getNumberOfContribs();
    if(nContributions==0) continue;
    int trackID = simHit.getContrib(0).trackID;  //There should only be 1 contributions for Hcal hits
    double charge = particleMap[trackID].getCharge();

    int absorber=0;  //no absorber at back section
    if(section!=ldmx::HcalID::HcalSection::BACK)
    {
      if(layer%2==1) absorber=2; //for odd layer numbers at side sections, absorber is at positive end   //TODO: Is this true?
      else
      {
        switch (section)
        {
          case ldmx::HcalID::HcalSection::TOP    : 
          case ldmx::HcalID::HcalSection::RIGHT  : absorber=2;  //absorber at positive end
                                                   break;
          case ldmx::HcalID::HcalSection::BOTTOM : 
          case ldmx::HcalID::HcalSection::LEFT   : absorber=-2;  //absorber at negative end
                                                   break;
        }
      }
    }

    photonGenerator->second->MakePhotons(pos1Local, pos2Local, simHit.getPreStepTime(), simHit.getPostStepTime(),
                                         simHit.getVelocity()/CLHEP::c_light, charge,
                                         simHit.getEdep(),
                                         simHit.getPathLength(),
                                         absorber);


    const std::vector<double> &timesPosEnd=photonGenerator->second->GetArrivalTimes(1);  //SiPM1 is at the positive and in the photonGenerator
    const std::vector<double> &timesNegEnd=photonGenerator->second->GetArrivalTimes(0);  //SiPM0 is at the negative end in the photonGenerator
    switch (section)
    {
      case ldmx::HcalID::HcalSection::BACK   : photonTimes[posEndID].insert(photonTimes[posEndID].end(),timesPosEnd.begin(),timesPosEnd.end());
                                               photonTimes[negEndID].insert(photonTimes[negEndID].end(),timesNegEnd.begin(),timesNegEnd.end());
                                               break;
      case ldmx::HcalID::HcalSection::TOP    : 
      case ldmx::HcalID::HcalSection::RIGHT  : photonTimes[negEndID].insert(photonTimes[negEndID].end(),timesNegEnd.begin(),timesNegEnd.end());
                                               break;
      case ldmx::HcalID::HcalSection::BOTTOM : 
      case ldmx::HcalID::HcalSection::LEFT   : photonTimes[posEndID].insert(photonTimes[posEndID].end(),timesPosEnd.begin(),timesPosEnd.end());
                                               break;
    }

#else

    //******************//
    //***  old code  ***//
    //******************//
    auto idh = hitsByID.find(hitID);
    if (idh == hitsByID.end()) {
      hitsByID[hitID] = std::vector<const ldmx::SimCalorimeterHit*>(1, &simHit);
    } else {
      idh->second.push_back(&simHit);
    }

#endif

  }

  /******************************************************************************************
   * HGCROC Emulation on Simulated Hits (grouped by HcalID)
   ******************************************************************************************/

#ifdef NEWCODE
  //loop over all SiPMs in order to also catch SiPMs that don't have hits, but may get dark counts
  for (int section = 0; section < hcalGeometry.getNumSections(); ++section)
  {
    for (int layer = 1; layer <= hcalGeometry.getNumLayers(section); ++layer)
    {
      for(int strip = 0; strip < hcalGeometry.getNumStrips(section, layer); ++strip)
      {
        typedef std::vector<ldmx::HgcrocDigiCollection::Sample> digisType;
        std::vector<std::pair<ldmx::HcalDigiID,digisType> > digisToAdd;   //can hold digis from both ends
        for(int end = 0; end <2; ++end)
        {
          ldmx::HcalDigiID channelID(section, layer, strip, end);

          if(section!=ldmx::HcalID::HcalSection::BACK)
          {
            if(layer%2==1)
            {
              if(end==0) continue; //for odd layer numbers at side sections, absorber is at pos end   //TODO: Is this true?
            }
            else
            {
              if(end==0 && section==ldmx::HcalID::HcalSection::TOP) continue;    //absorber at pos end
              if(end==0 && section==ldmx::HcalID::HcalSection::RIGHT) continue;  //absorber at pos end
              if(end==1 && section==ldmx::HcalID::HcalSection::BOTTOM) continue; //absorber at neg end
              if(end==1 && section==ldmx::HcalID::HcalSection::LEFT) continue;   //absorber at neg end
            }
          }

          if(randFlat_.fire() < deadSiPMProbability_) continue;  //assume that this random SiPM is dead
                                                                 //TODO: it may be better to select dead SiPMs at the beginning ot the run,
                                                                 //and keep these dead SiPMs for the entire run

          std::vector<std::pair<double,size_t> > photonTimesWithIndex;  //pair of photon time and index in the original photon vector
                                                                        //this is needed by the charge generator as used in the Mu2e CRV

          //find this HcalDigiID of this SiPM in the photonTimes map
          //and fill the new photon times/index vector
          auto times = photonTimes.find(channelID);
          if(times!=photonTimes.end())
          {
            for(size_t index=0; index<times->second.size(); ++index)
            {
              double time = times->second.at(index);
              photonTimesWithIndex.emplace_back(time,index);
            }
          }

          //generate individual pixes chargess based on arriving photons and add noise
          //return a charge vector of a type used by Mu2e
          std::vector<SiPMresponse> charges;
          chargeGenerator_->Simulate(photonTimesWithIndex, charges, digitizationStart_, digitizationEnd_);

          //convert vector of charges to vector of voltages
          std::vector<std::pair<double,double>> voltages;
          for(size_t index=0; index<charges.size(); ++index)
          {
            double voltage = charges.at(index)._chargeInPEs * singlePixelPeakVoltage_;
            double time = charges.at(index)._time;
            voltages.emplace_back(voltage,time);
          }

          //create digis //TODO: try to reduce having these digis copied twice
          digisType digis; 
          if(hgcroc_->digitize(channelID.raw(), voltages, digis)) digisToAdd.emplace_back(channelID,digis);  //can hold digis from both ends
        } //ends

        if(section == ldmx::HcalID::HcalSection::BACK)
        {
          if(digisToAdd.size()!=2) continue;  //need digis on both ends for back section
          hcalDigis.addDigi(digisToAdd[0].first.raw(),digisToAdd[0].second);
          hcalDigis.addDigi(digisToAdd[1].first.raw(),digisToAdd[1].second);
for(size_t i=0; i<digisToAdd[0].second.size(); ++i) std::cout<<"BACK 0 "<<i<<" "<<digisToAdd[0].second.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
for(size_t i=0; i<digisToAdd[1].second.size(); ++i) std::cout<<"BACK 1 "<<i<<" "<<digisToAdd[1].second.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
        }
        else
        {
          if(digisToAdd.size()!=1) continue;  //need digis on only one end for all other sections
          hcalDigis.addDigi(digisToAdd[0].first.raw(),digisToAdd[0].second);
for(size_t i=0; i<digisToAdd[0].second.size(); ++i) std::cout<<section<<"  "<<i<<" "<<digisToAdd[0].second.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
        }

      } //strips
    } //layers
  } //sections
std::cout<<"--------------------------------------------------------"<<std::endl;

#else
  for (auto const& simBar : hitsByID) {
    ldmx::HcalID detID(simBar.first);
    int section = detID.section();
    int layer = detID.layer();
    int strip = detID.strip();

    // get position
    double half_total_width = hcalGeometry.getHalfTotalWidth(section, layer);
    double ecal_dx = hcalGeometry.getEcalDx();
    double ecal_dy = hcalGeometry.getEcalDy();

    // contributions
    std::vector<std::pair<double, double>> pulses_posend;
    std::vector<std::pair<double, double>> pulses_negend;

    for (auto psimHit : simBar.second) {
      const ldmx::SimCalorimeterHit& simHit = *psimHit;

      std::vector<float> position = simHit.getPosition();

      /**
       * Define two pulses: with positive and negative ends.
       * For this we need to:
       * (1) Find the position along the bar:
       *     For back Hcal: x (y) for horizontal (vertical) layers.
       *     For side Hcal: x (top,bottom) and y (left,right).
       *
       * (2) Define the end of the bar:
       *     The end of an HcalDigiID is based on its distance (x,y) along the
       *     bar.
       *     - A positive end (endID=0), corresponds to top,left.
       *     - A negative end (endID=1), corresponds to bottom,right.
       *     For back Hcal:
       *     - if the position along the bar > 0, the close pulse's end is 0,
       *     else 1.
       *     For side Hcal:
       *     - if the position along the bar > half_width point of the bar, the
       *     close pulse's end is 0, else 1.
       *     The far pulse's end will be opposite to the close pulse's end.
       *
       * (3) Find the distance to each end (positive and negative) from the
       *     origin.
       *     For the back Hcal, the half point of the bar coincides with the
       *     coordinates of the origin.
       *     For the side Hcal, the length of the bar from the origin is:
       *     - 2 *(half_width) - Ecal_dx(y)/2 away from the positive end, and,
       *     - Ecal_dx(y) away from the negative end.
       */
      float distance_along_bar, distance_ecal;
      float distance_close, distance_far;
      int end_close;
      if (section == ldmx::HcalID::HcalSection::BACK) {
        distance_along_bar =
            hcalGeometry.layerIsHorizontal(layer) ? position[0] : position[1];
        end_close = (distance_along_bar > 0) ? 0 : 1;
        distance_close = half_total_width;
        distance_far = half_total_width;
      } else {
        if ((section == ldmx::HcalID::HcalSection::TOP) ||
            ((section == ldmx::HcalID::HcalSection::BOTTOM))) {
          distance_along_bar = position[0];
          distance_ecal = ecal_dx;
        } else if ((section == ldmx::HcalID::HcalSection::LEFT) ||
                   (section == ldmx::HcalID::HcalSection::RIGHT)) {
          distance_along_bar = position[1];
          distance_ecal = ecal_dy;
        }
        end_close = (distance_along_bar > half_total_width) ? 0 : 1;
        distance_close = (end_close == 0)
                             ? 2 * half_total_width - distance_ecal / 2
                             : distance_ecal / 2;
        distance_far = (end_close == 0)
                           ? distance_ecal / 2
                           : 2 * half_total_width - distance_ecal / 2;
      }

      // Calculate voltage attenuation and time shift for the close and far
      // pulse.
      float v = 299.792 /
                1.6;  // velocity of light in Polystyrene, n = 1.6 = c/v mm/ns
      double att_close =
          exp(-1. * ((distance_close - fabs(distance_along_bar)) / 1000.) /
              attlength_);
      double att_far =
          exp(-1. * ((distance_far + fabs(distance_along_bar)) / 1000.) /
              attlength_);
      double shift_close =
          fabs((distance_close - fabs(distance_along_bar)) / v);
      double shift_far = fabs((distance_far + fabs(distance_along_bar)) / v);

      // Get voltages and times.
      for (int iContrib = 0; iContrib < simHit.getNumberOfContribs();
           iContrib++) {
        double voltage = simHit.getContrib(iContrib).edep * MeV_;
        double time =
            simHit.getContrib(iContrib).time;  // global time (t=0ns at target)
        time -= position.at(2) /
                299.702547;  // shift light-speed particle traveling along z

        if (end_close == 0) {
          pulses_posend.emplace_back(voltage * att_close, time + shift_close);
          pulses_negend.emplace_back(voltage * att_far, time + shift_far);
        } else {
          pulses_posend.emplace_back(voltage * att_far, time + shift_far);
          pulses_negend.emplace_back(voltage * att_close, time + shift_close);
        }
      }
    }

    /**
     * Now we have all the sub-hits from all the simhits
     * Digitize:
     * For back Hcal return two digis.
     * For side Hcal we choose which pulse to readout based on
     * the position of the hit and the sub-section.
     * For Top and Left we read the positive end digi.
     * For Bottom and Right we read the negative end digi.
     **/
    if (section == ldmx::HcalID::HcalSection::BACK) {
      std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAddPosend,
          digiToAddNegend;
      ldmx::HcalDigiID posendID(section, layer, strip, 0);
      ldmx::HcalDigiID negendID(section, layer, strip, 1);
      if (hgcroc_->digitize(posendID.raw(), pulses_posend, digiToAddPosend) &&
          hgcroc_->digitize(negendID.raw(), pulses_negend, digiToAddNegend)) {
        hcalDigis.addDigi(posendID.raw(), digiToAddPosend);
        hcalDigis.addDigi(negendID.raw(), digiToAddNegend);
for(size_t i=0; i<digiToAddPosend.size(); ++i) std::cout<<"BACK 0"<<i<<" "<<digiToAddPosend.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
for(size_t i=0; i<digiToAddNegend.size(); ++i) std::cout<<"BACK 1"<<i<<" "<<digiToAddNegend.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
      }  // Back Hcal needs to digitize both pulses or none
    } else {
      bool is_posend = false;
      std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAdd;
      if ((section == ldmx::HcalID::HcalSection::TOP) ||
          (section == ldmx::HcalID::HcalSection::LEFT)) {
        is_posend = true;
      } else if ((section == ldmx::HcalID::HcalSection::BOTTOM) ||
                 (section == ldmx::HcalID::HcalSection::RIGHT)) {
        is_posend = false;
      }
      if (is_posend) {
        ldmx::HcalDigiID digiID(section, layer, strip, 0);
        if (hgcroc_->digitize(digiID.raw(), pulses_posend, digiToAdd)) {
          hcalDigis.addDigi(digiID.raw(), digiToAdd);
for(size_t i=0; i<digiToAdd.size(); ++i) std::cout<<section<<"  "<<i<<" "<<digiToAdd.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
        }
      } else {
        ldmx::HcalDigiID digiID(section, layer, strip, 1);
        if (hgcroc_->digitize(digiID.raw(), pulses_negend, digiToAdd)) {
          hcalDigis.addDigi(digiID.raw(), digiToAdd);
for(size_t i=0; i<digiToAdd.size(); ++i) std::cout<<section<<"  "<<i<<" "<<digiToAdd.at(i).adc_t()<<std::endl;
std::cout<<std::endl;
        }
      }
    }
  }
std::cout<<"--------------------------------------------------------"<<std::endl;

  /******************************************************************************************
   * Noise Simulation on Empty Channels
   *****************************************************************************************/
  if (noise_) {
    int numChannels = 0;
    for (int section = 0; section < hcalGeometry.getNumSections(); section++) {
      int numChannelsInSection = 0;
      for (int layer = 1; layer <= hcalGeometry.getNumLayers(section);
           layer++) {
        numChannelsInSection += hcalGeometry.getNumStrips(section, layer);
      }
      // for back Hcal we have double readout, therefore we multiply the number
      // of channels by 2.
      if (section == ldmx::HcalID::HcalSection::BACK) {
        numChannelsInSection *= 2;
      }
      numChannels += numChannelsInSection;
    }
    int numEmptyChannels = numChannels - hcalDigis.getNumDigis();
    // noise generator gives us a list of noise amplitudes [mV] that randomly
    // populate the empty channels and are above the readout threshold
    auto noiseHitAmplitudes{
        noiseGenerator_->generateNoiseHits(numEmptyChannels)};
    std::vector<std::pair<double, double>> fake_pulse(1, {0., 0.});
    for (double noiseHit : noiseHitAmplitudes) {
      // generate detector ID for noise hit
      // making sure that it is in an empty channel
      unsigned int noiseID;
      int sectionID, layerID, stripID, endID;
      do {
        sectionID = noiseInjector_->Integer(hcalGeometry.getNumSections());
        layerID = noiseInjector_->Integer(hcalGeometry.getNumLayers(sectionID));
        // set layer to 1 if the generator says it is 0 (geometry map starts
        // from 1)
        if (layerID == 0) layerID = 1;
        stripID = noiseInjector_->Integer(
            hcalGeometry.getNumStrips(sectionID, layerID));
        endID = noiseInjector_->Integer(2);
        if ((sectionID == ldmx::HcalID::HcalSection::TOP) ||
            (sectionID == ldmx::HcalID::HcalSection::LEFT)) {
          endID = 0;
        } else if ((sectionID == ldmx::HcalID::HcalSection::BOTTOM) ||
                   (sectionID == ldmx::HcalID::HcalSection::RIGHT)) {
          endID = 1;
        }
        auto detID = ldmx::HcalDigiID(sectionID, layerID, stripID, endID);
        noiseID = detID.raw();
      } while (hitsByID.find(noiseID) != hitsByID.end());
      hitsByID[noiseID] =
          std::vector<const ldmx::SimCalorimeterHit*>();  // mark this as used

      // get a time for this noise hit
      fake_pulse[0].second = noiseInjector_->Uniform(clockCycle_);

      // noise generator gives the amplitude above the readout threshold
      // we need to convert it to the amplitude above the pedestal
      double gain = hgcroc_->gain(noiseID);
      fake_pulse[0].first = noiseHit +
                            gain * hgcroc_->readoutThreshold(noiseID) -
                            gain * hgcroc_->pedestal(noiseID);

      if (sectionID == ldmx::HcalID::HcalSection::BACK) {
        std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAddPosend,
            digiToAddNegend;
        ldmx::HcalDigiID posendID(sectionID, layerID, stripID, 0);
        ldmx::HcalDigiID negendID(sectionID, layerID, stripID, 1);
        if (hgcroc_->digitize(posendID.raw(), fake_pulse, digiToAddPosend) &&
            hgcroc_->digitize(negendID.raw(), fake_pulse, digiToAddNegend)) {
          hcalDigis.addDigi(posendID.raw(), digiToAddPosend);
          hcalDigis.addDigi(negendID.raw(), digiToAddNegend);
        }
      } else {
        std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAdd;
        if (hgcroc_->digitize(noiseID, fake_pulse, digiToAdd)) {
          hcalDigis.addDigi(noiseID, digiToAdd);
        }
      }
    }  // loop over noise amplitudes
  }    // if we should add noise

#endif

  event.add(digiCollName_, hcalDigis);

  return;
}  // produce

}  // namespace hcal

DECLARE_PRODUCER_NS(hcal, HcalDigiProducer);
