#ifndef HCAL_HCALDIGIPRODUCER_H_
#define HCAL_HCALDIGIPRODUCER_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <memory>  //for smart pointers
#include <set>     //for tracking used detector IDs

//----------//
//   LDMX   //
//----------//
#include "DetDescr/HcalDigiID.h"
#include "DetDescr/HcalGeometry.h"
#include "DetDescr/HcalID.h"
#include "Framework/EventProcessor.h"
#include "Hcal/HcalPhotonGenerator.h"
#include "Hcal/HcalChargeGenerator.h"
#include "Recon/Event/EventConstants.h"
#include "Recon/Event/HgcrocDigiCollection.h"
#include "SimCore/Event/SimCalorimeterHit.h"
#include "Tools/HgcrocEmulator.h"
#include "Tools/NoiseGenerator.h"

namespace hcal {

/**
 * @class HcalDigiProducer
 * @brief Performs basic HCal digitization
 */
class HcalDigiProducer : public framework::Producer {
 public:
  /**
   * Constructor
   * Makes unique noise generator and injector for this class
   */
  HcalDigiProducer(const std::string& name, framework::Process& process);

  /// Default destructor
  virtual ~HcalDigiProducer() = default;

  /**
   * Configure this producer from the python configuration.
   * Sets event constants and configures the noise generator, noise injector,
   * and pulse function. Creates digi collection
   */
  void configure(framework::config::Parameters&) override;

  /**
   * Simulates measurement of pulse and creates digi collection for input event.
   */
  void produce(framework::Event& event) override;

 private:
  // Python Configuration Parameters

  /// input hit collection name
  std::string inputCollName_;

  /// input pass name
  std::string inputPassName_;

  /// output hit collection name
  std::string digiCollName_;

  /// Time interval for chip clock in ns
  double clockCycle_;

  /// Depth of ADC buffer.
  int nADCs_;

  /// Index for the Sample Of Interest in the list of digi samples
  int iSOI_;

  /// Conversion from energy in MeV to voltage in mV
  double MeV_;

  /// Strip attenuation length [m]
  double attlength_;

  // random number generators
  CLHEP::HepJamesRandom engine_;
  CLHEP::RandFlat       randFlat_;
  CLHEP::RandGaussQ     randGaussQ_;
  CLHEP::RandPoissonQ   randPoissonQ_;

  /****************************************************************************************** 
   * Usage of photon and charge generator 
   *
   * Photon generator: 
   *  generates individual photons based on the deposited energy of the track going through the scintillator
   *  uses lookup tables (as function of scintillator length) to determine Probability 
   * Charge generator: 
   *  simulates the response of the SiPM pixels to incoming photons
   *  generate individual pixes chargess based on arriving photons and add noise
   ******************************************************************************************/

  /// PHOTON GENERATOR
  /// mean scintillation yield
  double scintillationYield_;
  /// sigma of scintillation yield
  double scintillationYieldSigma_;
  /// cut-off for random scintillation yield
  double scintillationYieldCutoffLow_;
  double scintillationYieldCutoffHigh_;

  // map of photon generators per scintillatorLength and reflectorType
  std::map<std::pair<int,int>, std::shared_ptr<HcalPhotonGenerator> > photonGenerators_;

  // variables for charge generator
  double singlePixelPeakVoltage_;  // Peak voltage of the single pixel waveform [mV]
  double deadSiPMProbability_;

  std::shared_ptr<HcalChargeGenerator> chargeGenerator_;

  ///////////////////////////////////////////////////////////////////////////////////////
  // Other member variables

  /// Put noise into empty channels, not configurable, only helpful in
  /// development
  bool noise_{true};

  /// Use photon generator
  bool photongen_{true};

  /// Hgcroc Emulator to digitize analog voltage signals
  std::unique_ptr<ldmx::HgcrocEmulator> hgcroc_;

  /// Conversion from time in ns to ticks of the internal clock
  double ns_;

  /// Generates noise hits based off of number of cells that are not hit
  std::unique_ptr<ldmx::NoiseGenerator> noiseGenerator_;

  /// Generates Gaussian noise on top of real hits
  std::unique_ptr<TRandom3> noiseInjector_;
};
}  // namespace hcal

#endif
