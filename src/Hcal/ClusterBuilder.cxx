#include "Hcal/ClusterBuilder.h"

namespace hcal {

  ClusterGeometry::ClusterGeometry( ldmx::HcalGeometry hcalGeometry, int num_neighbors=4 ) {
    
    auto id_map = hcalGeometry.getStripPositionMap();
      
    // find adjacent strips in the same layer or adjacent layers
    for(auto pair1 = id_map.begin(); pair1 != id_map.end(); pair1++){
      for(auto pair2 = pair1; pair2 != id_map.end(); pair2++){
	if (pair1==pair2) continue;
	auto &id1 = pair1->first;
	auto &id2 = pair2->first;
	auto &pos1 = pair1->second;
	auto &pos2 = pair2->second;
	if (id1.section() != id2.section()) continue;
	if (id1.layer() == id2.layer()) {
	  if (fabs(id1.strip() - id2.strip()) == 1) {
	    AddStripNeighbor(id1.raw(), id2.raw());
	  }
	}
	if (fabs(id1.layer() - id2.layer()) == 1) {
	  // cout << "neigh layers, x,y,strip1 " <<  pos1.X() << " " << pos1.Y() << " " << id1.strip() << " x,y,strip2 " << pos2.X() << " " << pos2.Y() << " " << id2.strip() << endl;
	  // consider all neighboring strips in the map
	  float d = 0;
	  if (hcalGeometry.layerIsHorizontal(id1.layer())) {
	    // if a layer is horizontal, compute x distance
	    d = fabs(pos1.X() - pos2.X());
	  }
	  else {
	    d = fabs(pos1.Y() - pos2.Y());
	  }
	  // width of scintillator is 50mm
	  if ( d <= num_neighbors*50*2) {
	    AddLayerNeighbor(id1.raw(), id2.raw());
	  }
	}
	
      } // end loop over pair2
    } // end loop over pair1
    
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
      
      if (isLocalMax && (hit.e > seed_threshold_2d_) && !hit.used){
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
    while(i_neighbor < num_neighbors_){
      
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
	       hits_by_id[n].e > neighbor_threshold_2d_ ){
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
	  // double w = std::max(0., log(h.e/MIN_ENERGY_)); // use log-e wgt
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
    layer_clusters.resize(layer_max_); // first 100 layers
    for(auto &clus : all_clusters){
      layer_clusters[clus.layer].push_back(clus);
    }
    
    // sort clusters in each layer by energy
    for(auto &clusters : layer_clusters) 
      ESort(clusters);
    
    
    if(debug){
      cout << "--------\n3d: sorted 2d inputs" << endl;
      for(auto &clusters : layer_clusters) 
	for(auto &c : clusters) c.Print();
      cout << " size " << layer_clusters.size() << endl;
    }
    
    // identify cluster seed layer (max shower)
    // repeat, if not all clusters are used
    bool building=true;
    std::vector<Cluster> clusters3d;
    while(building){

      Cluster cluster3d;
      cluster3d.is2D = false;

      // find seed, from global shower max
      // derived from the first available cluster in each layer
	
      int layer_showermax = -1;
      double e_showermax = seed_threshold_3d_;
      for(auto &clusters : layer_clusters) {
	if(clusters.size() > 0){
	  if(clusters.at(0).e > e_showermax) {
	    layer_showermax = clusters.at(0).layer;
	    e_showermax = clusters.at(0).e;
	  }
	}
      }

      if (layer_showermax==-1)
	break;

      if(debug) 
	cout << "e_showermax " << e_showermax << " layer " << layer_showermax << endl;

      auto &clusters2d=layer_clusters[layer_showermax];
      if(debug){
	clusters2d[0].Print();
      }
      
      cluster3d.clusters2d={ clusters2d[0] };
      cluster3d.depth = 1;
      cluster3d.first_layer = layer_showermax;
      cluster3d.last_layer = layer_showermax;
      clusters2d.erase( clusters2d.begin() );

      // loop to front of hcal from the showermax, and then to back of hcal	
      for( int ilayer=1; ilayer<layer_max_; ilayer++) {
	int test_layer = ilayer;
	if (ilayer < layer_max_ - layer_showermax)
	  test_layer = layer_showermax+ilayer;
	else
	  test_layer = layer_max_-ilayer-1;

	if(layer_clusters[test_layer].size() <= 0) continue;
	
	auto &last_seed2d = cluster3d.clusters2d.back().seed;
	if(test_layer == layer_showermax-1) 
	  last_seed2d = cluster3d.clusters2d.front().seed;
	
	auto &clusters2d=layer_clusters[test_layer];

	if(debug) 
	  cout << "look at clusters in layer " << test_layer << endl;
	
	for(int iclus2d=0; iclus2d<clusters2d.size(); iclus2d++){
	  // require an energy threshold for the 2d clusters
	  if( clusters2d[iclus2d].e <= neighbor_threshold_3d_ )
	    continue;
	  
	  // check if 2d seed is neighboring
	  // decided on dx and dy distances if TOA info is not used
	  // or full dr = sqrt(dx^2+dy^2) if TOA info is used
	  auto &seed2d = clusters2d[iclus2d].seed;

	  if(!geom->CheckLayerNeighbor(last_seed2d,seed2d)) {
	    if(debug){
	      cout << "  -- " << iclus2d << " seed " << seed2d << endl; 
	      cout << " not neigh " << endl;
	      clusters2d[iclus2d].Print();
	    }
	  }
	  
	  if(last_seed2d==seed2d || 
	     geom->CheckLayerNeighbor(last_seed2d,seed2d)){
	    
	    if(debug){
	      cout << " extend: ";
	      clusters2d[iclus2d].Print();
	    }
	    if(test_layer < cluster3d.first_layer)
	      cluster3d.first_layer = test_layer;
	    if(test_layer > cluster3d.last_layer)
	      cluster3d.last_layer = test_layer;
			    
	    cluster3d.clusters2d.push_back(clusters2d[iclus2d]);
	    cluster3d.depth++;
	    clusters2d.erase( clusters2d.begin()+iclus2d );
	    break;			    
	  }
	}
      } // end loop over layers

      if(debug) {
	for(auto &clusters : layer_clusters)
	  for(auto &c : clusters) c.Print();
      }
      
      if(cluster3d.depth==0){
	building=false;
      }
      else {
	if(cluster3d.depth >= 2)
	  clusters3d.push_back(cluster3d);
      }

    } // building
  
    // post-process
    for(auto &c : clusters3d){
      c.e=0;
      c.x=0;
      c.y=0;
      c.z=0;
      c.xx=0;
      c.yy=0;
      c.zz=0;
      float sumw=0;
      c.hits.clear();
      for(auto c2 : c.clusters2d){
	c.e += c2.e;
	//cout << "3d: " << c2.e << " " << log(c2.e/MIN_TP_ENERGY) << endl;
	double energy = c2.e;
	double w = std::max(0., energy);
	c.x += c2.x * w;
	c.y += c2.y * w;
	c.z += c2.z * w;
	c.xx += c2.x * c2.x * w;
	c.yy += c2.y * c2.y * w;
	c.zz += c2.z * c2.z * w;
	sumw += w;
	for(const auto &hit : c2.hits){
	  c.hits.push_back(hit);
	}
      }
      c.x /= sumw;
      // cout << "x: " << c.x << endl;
      c.y /= sumw;
      c.z /= sumw;
      c.xx /= sumw; //now is <x^2>
      c.yy /= sumw;
      c.zz /= sumw;
      c.xx = sqrt(c.xx - c.x * c.x); //now is sqrt(<x^2>-<x>^2)
      c.yy = sqrt(c.yy - c.y * c.y); 
      c.zz = sqrt(c.zz - c.z * c.z);
    }
    
    if(debug){
      cout << "--------\nFound 3d clusters" << endl;
      for(auto &c : clusters3d) c.Print3d();
    }

    all_clusters.clear();
    all_clusters.insert(all_clusters.begin(), 
			clusters3d.begin(), clusters3d.end());
  }
  
  void ClusterBuilder::BuildClusters(){
    if(debug){
      cout << "--------\nAll hits" << endl;
      for(auto &hit : all_hits) hit.Print();
    }

    // cluster the hits in each layer
    Build2DClusters();

    // cluster the hits across layers
    Build3DClusters();

    // sort by energy
    ESort(all_clusters);
  }
    
}
