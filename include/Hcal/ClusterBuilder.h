#ifndef CLUSTERBUILDER_H
#define CLUSTERBUILDER_H

#include <map>
#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "TFitResult.h"
#include "TVector3.h"

#include "DetDescr/HcalID.h"
#include "DetDescr/HcalGeometry.h"
#include "Hcal/Event/HcalCluster.h"

using std::cout;
using std::endl;

namespace hcal {

  class ClusterGeometry{
  public:
    std::map<int, std::vector<int> > strip_neighbors; 
    std::map<int, std::vector<int> > layer_neighbors; 

    ClusterGeometry( ldmx::HcalGeometry hcalGeometry, int num_neighbors);
      
    // add a strip neighbor
    void AddStripNeighbor(int id1, int id2) {
      if( strip_neighbors.count(id1) ) strip_neighbors[id1].push_back(id2);
      else strip_neighbors[id1]={id2};
      if( strip_neighbors.count(id2) ) strip_neighbors[id2].push_back(id1);
      else strip_neighbors[id2]={id1};
    };
    // are these neigboring strips?
    bool CheckStripNeighbor(int id1, int id2) {
      auto &ns = strip_neighbors[id1];
      return std::find(ns.begin(), ns.end(), id2) != ns.end();
    };
    // add a layer neighbor
    void AddLayerNeighbor(int id1, int id2) {
      if( layer_neighbors.count(id1) ) layer_neighbors[id1].push_back(id2);
      else layer_neighbors[id1]={id2};
      if( layer_neighbors.count(id2) ) layer_neighbors[id2].push_back(id1);
      else layer_neighbors[id2]={id1};
    };
    // are these neigboring layers?
    bool CheckLayerNeighbor(int id1, int id2) {
      auto &ns = layer_neighbors[id1];
      return std::find(ns.begin(), ns.end(), id2) != ns.end();
    };
    
  };

  // @TODO: Switch from Hit and Cluster classes to ldmx::HcalHit and ldmx::HcalCluster
  class Hit {
  public:
    float x=0;
    float y=0;
    float z=0;
    float e=0;
    ldmx::HcalID id;
    int rawid=-1;
    int layer=-1;
    int section=-1;
    int strip=-1;
    bool used=false;

    void Print(){
      cout << "Hit ("
	   << "e= " << e
	   << ", rawid=" << rawid 
	   << ", section= " << section
	   << ", layer= " << layer
           << ", strip= " << strip
	   << ", x= " << x
	   << ", y= " << y
	   << ", z= " << z
	   << ", used= " << used
	   << ")" 
	   << endl;
    }
  };
  class Cluster {
  public:
    std::vector<Hit> hits;
    std::vector<int> strips;
    std::vector<Cluster> clusters2d; // for 3d

    float x=0;
    float y=0;
    float z=0;
    float xx=0; // x RMS
    float yy=0; // y RMS
    float zz=0; // z RMS
    float e=0;
    int seed=-1; // raw id of seed
    int layer=-1;
    int section=-1;
    
    bool is2D=true;
    int depth=0;
    int first_layer=-1;
    int last_layer=-1;
    
    void Print() {
      cout << "Cluster ("
	   << "e= " << e
	   << ", seed id=" << seed
	   << ", x= " << x
	   << ", y= " << y
	   << ", z= " << z
	   << ", layer= " << layer
	   << ", nHit= " << hits.size()
	   << ", avgStrip= " << accumulate(strips.begin(), strips.end(), 0.0) / strips.size()
	   << ", seedStrip= " << strips.at(0)
	   << ")" 
	   << endl;
    }

    void Print3d(){
      cout << "Cluster ("
	   << "e= " << e
	   << ", seed id=" << seed
	   << ", x= " << x
	   << ", y= " << y
	   << ", z= " << z
	   << ", n2dClus= " << clusters2d.size()
	   << ", first_layer=" << first_layer
	   << ", depth=" << depth
	   << ")" 
	   << endl;
    }
    
    void PrintHits(){
      Print();
      for(auto &h : hits){
	cout << "  "; 
	h.Print();
      }
    }
    
  };
    
  class ClusterBuilder {
  public:
    std::vector<Hit> all_hits{};
    std::vector<Cluster> all_clusters{};
    std::vector<Cluster> clusters_2d{};
    
    // cluster geometry
    ClusterGeometry* geom;

    // debug
    bool debug = false;

    // energy thresholds (MeV)
    double seed_threshold_2d_;
    double neighbor_threshold_2d_;
    double seed_threshold_3d_;
    double neighbor_threshold_3d_;
    
    // number of neighboring strips
    int num_neighbors_;

    // maximum number of layers for a 3d cluster
    int layer_max_ = 100;

    // min energy
    double MIN_ENERGY_ = 0.5; // MeV

    // Set thresholds and number of neighbors
    void SetThresholds( double seed_threshold, double neighbor_threshold) {
      seed_threshold_ = seed_threshold;
      neighbor_threshold_ = neighbor_threshold;
    }
    void SetNeighbors ( int num_neighbors ) {
      num_neighbors_ = num_neighbors;
    }


    // add hit to all_hits
    void AddHit(ldmx::HcalHit h){
      Hit hit;
      hit.rawid = h.getID();
      ldmx::HcalID id(hit.rawid);
      hit.id = id;
      hit.e = h.getEnergy();
      hit.x = h.getXPos();
      hit.y = h.getYPos();
      hit.z = h.getZPos();
      hit.layer = h.getLayer();
      hit.strip = h.getStrip();
      hit.section = h.getSection();
      all_hits.push_back(hit);
    }

    // get clusters
    std::vector<Cluster> Get2DClusters(){return clusters_2d;}
    std::vector<Cluster> Get3DClusters(){return all_clusters;}

    // set cluster geometry
    void SetClusterGeo(ClusterGeometry* _geom){geom=_geom;}

    virtual void BuildClusters();
    std::vector<Cluster> Build2DClustersPerLayer( std::vector<Hit> hits);
    void Build2DClusters();
    void Build3DClusters();
  };

  template<class T>
    void ESort(std::vector<T> &v){
    std::sort( v.begin(),
		   v.end(),
	       [](const auto& lhs,const auto& rhs){
		 return lhs.e > rhs.e;
	       });
  }
  
} // namespace hcal

#endif
