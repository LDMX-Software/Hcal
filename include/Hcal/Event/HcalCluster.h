/**
 * @file HcalCluster.h
 * @brief Class that stores cluster information from the ECal
 * @author Josh Hiltbrand, University of Minnesota
 * @author Soophie Middleton, Caltech
 */

#ifndef EVENT_HCALCLUSTER_H_
#define EVENT_HCALCLUSTER_H_

// ROOT
#include "TObject.h"  //For ClassDef
#include "TString.h"

// STL
#include <iostream>
#include <set>

// ldmx-sw
#include "Hcal/Event/HcalHit.h"

namespace ldmx {

/**
 * @class HcalCluster
 * @brief Stores cluster information from the ECal
 */
class HcalCluster {
 public:
  /**
   * Class constructor.
   */
  HcalCluster();

  /**
   * Class destructor.
   */
  virtual ~HcalCluster();

  /**
   * Print a description of this object.
   */
  void Print() const;

  /**
   * Reset the HcalCluster object.
   */
  void Clear();

  /**
   * Take in the hits that make up the cluster.
   * @param hit The digi hit's entry number in the events digi
   * collection.
   */
  void addHits(const std::vector<const ldmx::HcalHit*> hitsVec);

  /**
   * Sets total energy for the cluster.
   * @param energy The total energy of the cluster.
   */
  void setEnergy(double energy) { energy_ = energy; }

  void setSeedEnergy(double energy) { seedenergy_ = energy;}
  
  /**
   * Sets total number of hits in the cluster.
   * @param nHits The total number of hits.
   */
  void setNHits(int nHits) { nHits_ = nHits; }

  /**
   * Sets total number of 2d clusters.
   * @param nClusters
   */
  void setN2DClusters(int n2d) { n2D_ = n2d; }
  
  /**
   * Adds strips to the cluster
   * @param stripsVec Vector of strips
   */
  void addStrips(const std::vector<int> stripsVec);

  void addStripsPerLayer(const std::vector<std::vector<int>> stripsVec);

  void setStripsOdd(const std::vector<int> stripsVec) {
    strips_odd_layer_ = stripsVec.size();
  }
  void setStripsEven(const std::vector<int> stripsVec) {
    strips_even_layer_ = stripsVec.size();
  }
  
  /**
   * Sets a sorted vector for the IDs of the hits
   * that make up the cluster.
   * @param IDs Sorted vector of hit IDs.
   */
  void setIDs(std::vector<unsigned int>& hitIDs) { hitIDs_ = hitIDs; }

  /**
   * Sets the three coordinates of the cluster centroid
   * @param x The x coordinate.
   * @param y The y coordinate.
   * @param z The z coordinate.
   */
  void setCentroidXYZ(double x, double y, double z) {
    centroidX_ = x;
    centroidY_ = y;
    centroidZ_ = z;
  }

  void setRMSXYZ(double xx, double yy, double zz) {
    rmsX_ = xx;
    rmsY_ = yy;
    rmsZ_ = zz;
  }

  void setLayer(int layer) { layer_ = layer; }
  
  void setTime(double x) { time_ = x; }

  void setDepth(double d) { depth_ = d; }
  
  /////////////////////////////////////////////

  // energy of cluster
  double getEnergy() const { return energy_; }

  // energy of cluster seed
  double getSeedEnergy() const { return seedenergy_; }
  
  // number of hits - equivalent to number of strips
  int getNHits() const { return nHits_; }

  // number of 2d clusters (for 3d clusters)
  int getN2DClusters() const { return n2D_; }

  // get layer (for 2d clusters)
  // or seedlayer (for 3d clusters)
  int getLayer() const { return layer_; }
  
  // position (weighted by energy)
  double getCentroidX() const { return centroidX_; }
  double getCentroidY() const { return centroidY_; }
  double getCentroidZ() const { return centroidZ_; }
  double getRMSX() const { return rmsX_; }
  double getRMSY() const { return rmsY_; }
  double getRMSZ() const { return rmsZ_; }

  // time (unused)
  double getTime() const { return time_; }

  // depth or number of layers in cluster (for 3d clusters)
  double getDepth() const { return depth_; }

  // get hit rawIDs (unused)
  const std::vector<unsigned int>& getHitIDs() const { return hitIDs_; }

  // get strips
  const std::vector<int>& getStrips() const { return strips_; }

  const std::vector< std::vector<int >> & getStripsPerLayer() const { return strips_per_layer_; }
  
  int getOddLayerStrips() const { return strips_odd_layer_; }

  int getEvenLayerStrips() const { return strips_even_layer_; }
  
  bool operator<(const HcalCluster& rhs) const {
    return this->getEnergy() < rhs.getEnergy();
  }

 private:
  std::vector<unsigned int> hitIDs_;
  std::vector<std::vector<int>> strips_per_layer_;
  std::vector<int> strips_;
  double energy_{0};
  double seedenergy_{0};
  int nHits_{0};
  int n2D_{0};
  int layer_{0};
  int strips_even_layer_{0};
  int strips_odd_layer_{0};
  double centroidX_{0};
  double centroidY_{0};
  double centroidZ_{0};
  double rmsX_{0};
  double rmsY_{0};
  double rmsZ_{0};
  double time_{0};
  double depth_{0};

  ClassDef(HcalCluster, 1);
};
}  // namespace ldmx

#endif
