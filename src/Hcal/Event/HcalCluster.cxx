#include "Hcal/Event/HcalCluster.h"

ClassImp(ldmx::HcalCluster)

    namespace ldmx {
  HcalCluster::HcalCluster() {}

  HcalCluster::~HcalCluster() { Clear(); }

  void HcalCluster::Print() const {
    std::cout << "HcalCluster { "
              << "Energy: " << energy_ << ", "
              << "Number of hits: " << nHits_ << " }" << std::endl;
  }

  void HcalCluster::Clear() {
    hitIDs_.clear();
    strips_.clear();
    strips_per_layer_.clear();
    energy_ = 0;
    nHits_ = 0;
    n2D_ = 0;
    layer_ = 0;
    depth_ = 0;
    time_ = 0;
    centroidX_ = 0;
    centroidY_ = 0;
    centroidZ_ = 0;
    rmsX_ = 0;
    rmsY_ = 0;
    rmsZ_ = 0;
    DXDZ_ = 0;
    DYDZ_ = 0;
    errDXDZ_ = 0;
    errDYDZ_ = 0;
  }

  void HcalCluster::addHits(const std::vector<const HcalHit *> hitsVec) {
    std::vector<unsigned int> vecIDs;
    for (unsigned int iHit = 0; iHit < hitsVec.size(); iHit++) {
      vecIDs.push_back(hitsVec[iHit]->getID());
    }
    setIDs(vecIDs);
  }

  void HcalCluster::addStrips(const std::vector<int> stripsVec) {
    // add only unique strips to a list
    strips_.clear();
    for (auto strip : stripsVec) {
      if (!std::count(strips_.begin(), strips_.end(), strip)) {
        strips_.push_back(strip);
      }
    }
  }

  void HcalCluster::addStripsPerLayer(const std::vector<std::vector<int>> stripsVec) {
    // strips_per_layer_.clear();
    std::vector<std::vector<int>> strips_per_layer(100);
    for (unsigned int layer=0; layer < stripsVec.size(); layer++) {
      for (auto strip : stripsVec.at(layer)) {
	strips_per_layer.at(layer).push_back(strip);
      }
    }
    strips_per_layer_ = strips_per_layer;
    // for (unsigned int l=0; l < strips_per_layer_.size(); l++) {
    //   for (auto strip : strips_per_layer_.at(l)) {
    // 	std::cout << " l " << l << " s " << strip << std::endl;
    //   }
    // }
    
  }
      

}  // namespace ldmx
