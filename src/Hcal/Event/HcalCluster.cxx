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
      
}  // namespace ldmx
