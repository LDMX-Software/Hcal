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

  /// energy threshold for 2D seed
  double seed_threshold_;
  /// energy threshold for 2D neighbors
  double neighbor_threshold_;
  /// number of 2D neighbors
  int num_neighbors_;
  
 private:

 public:
  HcalNewClusterProducer(const std::string& n, framework::Process& p)
      : Producer(n, p) {}

  virtual void configure(framework::config::Parameters& p) final override;
  virtual void produce(framework::Event& event) final override;

};  // HcalNewClusterProducer

}  // namespace hcal
#endif
