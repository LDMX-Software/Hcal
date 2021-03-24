#include "Hcal/Event/HcalCluster.h"

ClassImp(ldmx::HcalCluster)


std::ostream& operator<<(std::ostream& s, const ldmx::HcalCluster& hc) {
    return s << "HcalCluster { "
             << "Energy: " << hc.getEnergy() << ", "
             << "Position: (" << hc.getCentroidX() << ", "<< hc.getCentroidY() << ", "<< hc.getCentroidZ() << ") , "
             << "Number of hits: " << hc.getNHits() << " }" ;
  }
namespace ldmx {
  HcalCluster::HcalCluster() {}

  HcalCluster::~HcalCluster() { Clear(); }

  void HcalCluster::Print() const {
    std::cout << *this << std::endl;
  }

  void HcalCluster::Clear() {
    hitIDs_.clear();

    energy_ = 0;
    nHits_ = 0;
    centroidX_ = 0;
    centroidY_ = 0;
    centroidZ_ = 0;
  }

  void HcalCluster::addHits(const std::vector<const HcalHit *> hitsVec) {

    std::vector<unsigned int> vecIDs;
    for (int iHit = 0; iHit < hitsVec.size(); iHit++) {
      vecIDs.push_back(hitsVec[iHit]->getID());
    }

    setIDs(vecIDs);
  }
}  // namespace ldmx
