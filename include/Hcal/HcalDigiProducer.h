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
#include <TRandom2.h>

#include "DetDescr/HcalDigiID.h"
#include "DetDescr/HcalGeometry.h"
#include "DetDescr/HcalID.h"
#include "Framework/EventProcessor.h"
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

  enum class TimeSpreadType {
    GAUSSIAN = 0,
    UNIFORM = 1,
    CONSTANT = 2,
  };

  /** Create a randomized time delta based on the configured time spread type.
   * Currently defined options are:
   * 0: Gaussian -> parameters[0] = mean, parameters[1] = sigma
   * 1: Uniform -> parameters[0] = min, parameters[1] = max
   * 2: Constant -> parameters[0] = value
   */
  double makeTimeDelta(TimeSpreadType kind,
                       std::vector<double>& parameters) const;

  TimeSpreadType validateTimeSpreadType(int kind,
                                        std::vector<double>& parameters) const;

 private:
  ///////////////////////////////////////////////////////////////////////////////////////
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

  ///////////////////////////////////////////////////////////////////////////////////////
  // Other member variables

  /// Put noise into empty channels, not configurable, only helpful in
  /// development
  bool noise_{true};

  /// Hgcroc Emulator to digitize analog voltage signals
  std::unique_ptr<ldmx::HgcrocEmulator> hgcroc_;

  std::unique_ptr<TRandom2> randomTime_;
  bool doTimeSpreadPerSpill{false};
  TimeSpreadType timeSpreadPerSpillType{TimeSpreadType::UNIFORM};
  std::vector<double> timeSpreadPerSpillParameters{};

  bool doTimeSpreadPerHit{false};
  TimeSpreadType timeSpreadPerHitType{TimeSpreadType::GAUSSIAN};
  std::vector<double> timeSpreadPerHitParameters{};

  double flatTimeShift{0.};
  /// Conversion from time in ns to ticks of the internal clock
  double ns_;

  /// Generates noise hits based off of number of cells that are not hit
  std::unique_ptr<ldmx::NoiseGenerator> noiseGenerator_;

  /// Generates Gaussian noise on top of real hits
  std::unique_ptr<TRandom3> noiseInjector_;
};
}  // namespace hcal

#endif
