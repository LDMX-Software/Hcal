#ifndef HCAL_HCALNEWCLUSTERPRODUCER_H_
#define HCAL_HCALNEWCLUSTERPRODUCER_H_

#include "DetDescr/HcalGeometry.h"
#include "DetDescr/HcalID.h"

#include "Framework/EventDef.h"
#include "Framework/EventProcessor.h"

namespace hcal {

class HcalNewClusterProducer : public framework::Producer {
  /// name of pass of rechits to use
  std::string pass_name_{""};
  /// name of rechits to reconstruct
  std::string coll_name_{"HcalRecHits"};
  /// name of 2d clusters to reconstruct
  std::string cluster2d_coll_name_{"Hcal2DClusters"};
  /// name of 3d clusters to reconstruct
  std::string cluster3d_coll_name_{"Hcal3DClusters"};

  /// noise energy threshold for hits entering clustering
  double noise_threshold_;
  /// energy threshold for 2D/3D seed
  double seed_threshold_2d_;
  double seed_threshold_3d_;
  /// energy threshold for 2D/3D neighbors
  double neighbor_threshold_2d_;
  double neighbor_threshold_3d_;
  /// max xy
  double max_xy_2d_;
  double max_xy_3d_;
  /// number of 2D neighbors
  int num_neighbors_;
  /// whether to use TOA info or not
  bool use_toa_;
  
 private:

 public:
  HcalNewClusterProducer(const std::string& n, framework::Process& p)
      : Producer(n, p) {}

  virtual void configure(framework::config::Parameters& p) final override;
  virtual void produce(framework::Event& event) final override;

};  // HcalNewClusterProducer

}  // namespace hcal
#endif
