/**
 * @file HcalMIPTracking.h
 * @brief Class that performs MIP tracking in the Hcal
 * @author Matthew Solt, UVA
 */

#ifndef HCAL_HCALMIPTRACKING_H_
#define HCAL_HCALMIPTRACKING_H_

// ROOT
#include "TRandom3.h"
#include "TString.h"

// LDMX
#include "DetDescr/DetectorID.h"
#include "DetDescr/HcalID.h"
#include "Framework/Configure/Parameters.h"
#include "Framework/EventDef.h"
#include "Framework/EventProcessor.h"

namespace hcal {

/**
 * @class HcalMIPTracking
 * @brief Performs MIP tracking in the Hcal
 */
class HcalMIPTracking : public framework::Producer {
 public:
  HcalMIPTracking(const std::string& name, framework::Process& process);

  virtual ~HcalMIPTracking() { ; }

  /**
   * Configure the processor using the given user specified parameters.
   *
   * @param parameters Set of parameters used to configure this processor.
   */
  void configure(framework::config::Parameters& parameters) final override;

  virtual void produce(framework::Event& event);

 private:

  bool verbose_{false};

  double BACK_HCAL_START_Z_{840.};

  double MIP_MIN_PE_{36.};

  double MIP_MAX_PE_{900.};

  int MIN_TRACK_HITS_{4};

  int MIN_SEED_HITS_{3};

  double MAX_SEED_HIT_ERROR_{300.};

  double MAX_TRACK_EXTRAP_SIGMA_{5.};

  int MAX_LAYERS_CONSEC_MISSED_{2};

  double strip_position_resolution_{150.};

  int STRIPS_BACK_PER_LAYER_{60};

  int NUM_BACK_HCAL_LAYERS_{150};

  int STRIPS_SIDE_TB_PER_LAYER_{6};

  int NUM_SIDE_TB_HCAL_LAYERS_{31};

  int STRIPS_SIDE_LR_PER_LAYER_{31};

  int NUM_SIDE_LR_HCAL_LAYERS_{63};

  int SUPER_STRIP_SIZE_{1};

  std::vector<ldmx::HcalHit> FindIsolatedHits(std::vector<ldmx::HcalHit>& hits);

  float CalcDist(ldmx::HcalHit& h, std::vector<ldmx::HcalHit> &hits);

  std::vector<ldmx::HcalHit> chooseSeed(std::vector<std::vector<ldmx::HcalHit>> &seedlist);

  std::vector<std::vector<ldmx::HcalHit>> FindTracks(std::map<int, std::vector<ldmx::HcalHit>> &hitmap);

  std::vector<std::vector<ldmx::HcalHit>> FindSeeds(std::map<int, std::vector<ldmx::HcalHit>> &hitmap, std::vector<ldmx::HcalHit> &rmhit, int &index);

  std::vector<std::vector<ldmx::HcalHit>> CleanSeedList(std::vector<std::vector<ldmx::HcalHit>> &seedlist);

  bool isOnTrack(ldmx::HcalHit &hit, std::vector<ldmx::HcalHit> &rmhit);

  bool IsTriggered(std::map<int, std::vector<ldmx::HcalHit>> &hitmap);

  float vectorMean(std::vector<float>& vec);

  float vectorSD2(std::vector<float>& vec);

  float vectorSS(std::vector<float>& vec1, std::vector<float>& vec2);

  float* fitTrackLS(std::vector<ldmx::HcalHit> &hitlist);
};

}  // namespace hcal

#endif
