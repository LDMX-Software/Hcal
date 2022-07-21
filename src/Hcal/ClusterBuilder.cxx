#include "Hcal/ClusterBuilder.h" 

namespace hcal {

  ClusterGeometry::ClusterGeometry( std::map<ldmx::HcalID, TVector3> id_map ){
    // find adjacent strips in the same layer or adjacent layers
    for(auto pair1 = id_map.begin(); pair1 != id_map.end(); pair1++){
      for(auto pair2 = pair1; pair2 != id_map.end(); pair2++){
	if (pair1==pair2) continue;
	auto &id1 = pair1->first;
	auto &id2 = pair2->first;
	if (id1.section() != id2.section()) continue;
	if (id1.layer() == id2.layer()) {
	  if (fabs(id1.strip() - id2.strip()) == 1) {
	    AddStripNeighbor(id1.raw(), id2.raw());
	  }
	}
	if (fabs(id1.layer() - id2.layer()) == 1) {
	  AddLayerNeighbor(id1.raw(), id2.raw());
	}
      }
    }
    
  }
  
  std::vector<Cluster> ClusterBuilder::Build2DClustersPerLayer(std::vector<Hit> hits) {
    // map hits by id
    std::map<int,Hit> hits_by_id;
    for(auto &hit : hits) hits_by_id[hit.rawid]=hit;

    if(debug){
      cout << "--------\nBuild2DClustersLayer Input Hits" << endl;
      for(auto &hitpair : hits_by_id) hitpair.second.Print();
      cout << " ----- " << endl;
    }

    // find seeds
    std::vector<Cluster> clusters;
    for(auto &hitpair : hits_by_id){
      auto &hit = hitpair.second;
      bool isLocalMax = true;
      for(auto n : geom->strip_neighbors[hit.rawid]){
	// if(debug) {
	//   if (hits_by_id.count(n)) {
	//     cout << "neighbor is in list of hits " << hits_by_id[n].e << " rawid " << hits_by_id[n].rawid << " hit e " << hit.e << " raw id " << hit.rawid << endl;
        //  }
	// }
	if (hits_by_id.count(n) && hits_by_id[n].e > hit.e && hits_by_id[n].rawid != hit.rawid) {
	  isLocalMax = false;
	}
      }
      
      if (isLocalMax && (hit.e > seed_threshold) && !hit.used){
	hit.used=true;

	// if(debug) {
	//   cout << "local max " << endl;
	//   hit.Print();
	// }
	
	Cluster c;
	c.hits.push_back( hit );
	c.e = hit.e;
	c.x = hit.x;
	c.y = hit.y;
	c.z = hit.z;
	c.seed = hit.rawid;
	c.layer = hit.layer;
	c.section = hit.section;
	c.strips.push_back( hit.strip );
	clusters.push_back(c);
      }
    } // end loop over hits by ID
    
    if(debug){
      cout << "After seed finding: " << endl;
      for(auto &hitpair : hits_by_id) hitpair.second.Print();
      cout << "Print cluster" << endl;
      for(auto &c : clusters) c.Print();
      cout << " ----- " << endl;
    }

    // Add neighbors up to the specified limit
    int i_neighbor=0;
    while(i_neighbor < num_neighbors){
      
      // find unused neighbors for all clusters
      std::map<int, std::vector<int> > assoc_clus2hitIDs;
      int unused_hits = 0;
      for(int iclus=0; iclus<clusters.size(); iclus++){
	auto &clus=clusters[iclus];
	std::vector<int> neighbors;
	for(const auto &hit : clus.hits){
	  for(auto n : geom->strip_neighbors[hit.rawid]){
	    if(hits_by_id.count(n) && 
	       !hits_by_id[n].used && 
	       hits_by_id[n].e > neighbor_threshold ){
	      neighbors.push_back(n);
	      unused_hits ++;
	    }
	  }
	}
	assoc_clus2hitIDs[iclus] = neighbors;
      }

      if(unused_hits == 0) break;
      
      // check to how many clusters, each hit is associated to
      std::map<int, std::vector<int> > assoc_hitID2clusters;
      for(auto clus2hitID : assoc_clus2hitIDs){
	auto iclus = clus2hitID.first;
	auto &hitIDs = clus2hitID.second;
	for(const auto &hitID : hitIDs){
	  if( assoc_hitID2clusters.count(hitID) )
	    assoc_hitID2clusters[hitID].push_back(iclus);
	  else
	    assoc_hitID2clusters[hitID] = {iclus};
	}
      }

      // add associated hits to clusters
      for(auto hitID2clusters : assoc_hitID2clusters){
	auto hitID = hitID2clusters.first;
	auto iclusters = hitID2clusters.second;
	if (iclusters.size() == 0) continue;
	if (iclusters.size() == 1){
	  // if a hit is associated to only one cluster, simply add it to the cluster
	  auto &hit = hits_by_id[hitID];
	  auto iclus = iclusters[0];
	  hit.used=true;
	  clusters[ iclus ].hits.push_back( hit );
	  clusters[ iclus ].strips.push_back( hit.strip );
	  clusters[ iclus ].e += hit.e;
	} else {
	  // if many hits are associated to one cluster, do WINNER TAKES IT ALL
	  auto &hit = hits_by_id[hitID];
	  hit.used=true;
	  float maxE=0;
	  int maxE_idx=-1;
	  for (auto iclus : iclusters){
	    if(clusters[iclus].e > maxE){
	      maxE = clusters[iclus].e;
	      maxE_idx = iclus;
	    }
	  }
	  clusters[ maxE_idx ].hits.push_back(hit);
	  clusters[ maxE_idx ].strips.push_back( hit.strip );
	  clusters[ maxE_idx ].e += hit.e;
	}
      }

      // rebuild cluster properties based on the new hits
      // do log energy weighting here
      for(auto &c : clusters){
	c.e=0;
	c.x=0;
	c.y=0;
	c.z=0;
	c.xx=0;
	c.yy=0;
	c.zz=0;
	double sumw=0;
	for(auto h : c.hits){
	  auto hit=h;
	  c.e += h.e;
	  double energy = h.e;
	  // double w = std::max(0., log(h.e/MIN_ENERGY)); // use log-e wgt
	  double w = std::max(0., energy);
	  c.x += h.x * w;
	  c.y += h.y * w;
	  c.z += h.z * w;
	  c.xx += h.x * h.x * w;
	  c.yy += h.y * h.y * w;
	  c.zz += h.z * h.z * w;
	  sumw += w;
	}
	c.x /= sumw;
	c.y /= sumw;
	c.z /= sumw;
	c.xx /= sumw; //now is <x^2>
	c.yy /= sumw;
	c.zz /= sumw;
	c.xx = sqrt(c.xx - c.x * c.x); //now is sqrt(<x^2>-<x>^2)
	c.yy = sqrt(c.yy - c.y * c.y); 
	c.zz = sqrt(c.zz - c.z * c.z);
      }

      i_neighbor ++;

      if(debug){
	cout << "--------\nAfter " << i_neighbor 
	     << " neighbors" << endl;
	for(auto &hitpair : hits_by_id) hitpair.second.Print();
	cout << "Print cluster " << endl;
	for(auto &c : clusters) c.Print();
	cout << " ----- " << endl;
      }

      
      
    } // neighbor loop
    
    return clusters;
  }
  
  void ClusterBuilder::Build2DClusters() {      
    // @TODO: for now only use back hcal
    // loop over sections
    for (int section = 0; section < 1; section++) {
    
      // map hits by layer
      std::map<int, std::vector<Hit> > layer_hits;
      for(const auto hit : all_hits){
	auto l = hit.layer;
	if (hit.section != section) continue;
	if (layer_hits.count(l)){
	  layer_hits[l].push_back(hit);
	} else {
	  layer_hits[l]={hit};
	}
      }

      // run clustering in each layer
      for(auto &pair : layer_hits){
	if(debug){
	  cout << "Found " << pair.second.size()
	       << " hits in layer " << pair.first << endl;
	}
	auto clus = Build2DClustersPerLayer(pair.second);
	all_clusters.insert(all_clusters.end(), 
			    clus.begin(), clus.end());
      }
    }

    // sort by energy
    clusters_2d	= all_clusters;
    ESort(clusters_2d);
  }

  void ClusterBuilder::Build3DClusters() {
    if(debug){
      cout << "--------\nBuilding 3d clusters" << endl;
    }

    // sort 2d clusters by layer
    std::vector<std::vector<Cluster> > layer_clusters;
    //layer_clusters.resize(layer_max); // first 50 layers
    for(auto &clus : all_clusters){
      layer_clusters[clus.layer].push_back(clus);
    }
    // sort by energy
    for(auto &clusters : layer_clusters) 
      ESort(clusters);

    
    if(debug){
      cout << "--------\n3d: sorted 2d inputs" << endl;
      for(auto &clusters : layer_clusters) 
	for(auto &c : clusters) c.Print();
    }
    
    // @TODO: pass clusters from layer 0 to end, starting with highest energy
    bool building=true;
    std::vector<Cluster> clusters3d;

    if(debug){
      cout << "--------\nFound 3d clusters" << endl;
      // for(auto &c : clusters3d) c.Print3d();
    }
    
  }
  
  void ClusterBuilder::BuildClusters(){
    if(debug){
      cout << "--------\nAll hits" << endl;
      for(auto &hit : all_hits) hit.Print();
    }

    // cluster the hits in each layer
    Build2DClusters();

    // cluster the hits across layers
    // Build3DClusters();

    // sort by energy
    ESort(all_clusters);
  }
    
}
