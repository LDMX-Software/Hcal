#include "Hcal/ClusterBuilder.h"

namespace hcal {

ClusterGeometry::ClusterGeometry(ldmx::HcalGeometry hcalGeometry) {
  /**
   * Build map of posible strip neighbors in the same layer
   * or possible layer neighbors, given raw ID.
   */

  auto id_map = hcalGeometry.getStripPositionMap();

  for (auto pair1 = id_map.begin(); pair1 != id_map.end(); pair1++) {
    for (auto pair2 = pair1; pair2 != id_map.end(); pair2++) {
      if (pair1 == pair2) continue;
      auto &id1 = pair1->first;
      auto &id2 = pair2->first;
      auto &pos1 = pair1->second;
      auto &pos2 = pair2->second;
      // select only IDs with same section
      // TODO: add neighbors from different sections
      if (id1.section() != id2.section()) continue;

      // find adjacent strips in the same layer
      if (id1.layer() == id2.layer()) {
        if (fabs(id1.strip() - id2.strip()) == 1) {
          AddStripNeighbor(id1.raw(), id2.raw());
        }
      }

      // find adjacent strips in adjacent layers
      if (fabs(id1.layer() - id2.layer()) == 1) {
        // from the geometry, we have no way of knowing which strips are
        // neighboring so, we take all of the strips as neighbors and take the
        // distance between two hits in the 3d clustering
        AddLayerNeighbor(id1.raw(), id2.raw());
      }

    }  // end loop over pair2
  }    // end loop over pair1
}

void ClusterBuilder::ReBuild2DCluster(std::vector<Cluster> & clusters) {
  for (auto &c : clusters) {
    c.e = 0;
    c.x = 0;
    c.y = 0;
    c.z = 0;
    c.xx = 0;
    c.yy = 0;
    c.zz = 0;
    double sumw = 0;
    for (auto h : c.hits) {
      auto hit = h;
      c.e += h.e;
      double energy = h.e;
      //double w = 1;                                                                                                                                                                                       
      double w = std::max(0., energy);
      if (energy_weight_ == 1) {
        double w = std::max(0., log(c.e / MIN_ENERGY_));
      }
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
    c.xx /= sumw;  // now is <x^2>                                                                                                                                                                          
    c.yy /= sumw;
    c.zz /= sumw;
    c.xx = sqrt(c.xx - c.x * c.x);  // now is sqrt(<x^2>-<x>^2)                                                                                                                                             
    c.yy = sqrt(c.yy - c.y * c.y);
    c.zz = sqrt(c.zz - c.z * c.z);
  }
}

void ClusterBuilder::ReBuild3DCluster(std::vector<Cluster> & clusters3d) {
  if(debug)
    cout << "re-building 3d clusters " << clusters3d.size() << endl;
  
  for (auto &c : clusters3d) {
    c.e = 0;
    c.x = 0;
    c.y = 0;
    c.z = 0;
    c.xx = 0;
    c.yy = 0;
    c.zz = 0;
    double sumw = 0;
    double sumw_x = 0;
    double sumw_y = 0;
    c.hits.clear();
    c.avg_layer = 0;
    c.strips_oddlayer.clear();
    c.strips_evenlayer.clear();

    std::vector<std::vector<int>> s_perlayer(100);

    for (auto c2 : c.clusters2d) {
      c.e += c2.e;
      //if (debug) {
      //  c2.Print();
      //}
      double energy = c2.e;
      // double w = 1;                                                                                                                                                                                                                                                                                                                                               
      double w = std::max(0., energy);
      if (energy_weight_ == 1) {
        double w = std::max(0., log(c2.e / MIN_ENERGY_));
      }
      if (use_toa_) {
        // if using TOA information
        // only estimate the x position from vertical layers (even for v13!)
        // and only estimate the y position from horizontal layers (odd for v13!)
        // NOTE: Change to depend on geometry parity..
        //if(c2.layer % 2 == 0) {                                                                                                                                                                                                                                                                                                                                    
        c.x += c2.x * w;
        c.xx += c2.x * c2.x * w;
        sumw_x+=w;
        //}                                                                                                                                                                                                                                                                                                                                                          
        //else {                                                                                                                                                                                                                                                                                                                                                     
        c.y += c2.y * w;
        c.yy += c2.y * c2.y * w;
        sumw_y+=w;
        //}                                                                                                                                                                                                                                                                                                                                                          
      }
      else {
        c.x += c2.x * w;
        c.y += c2.y * w;
        c.xx += c2.x * c2.x * w;
        c.yy += c2.y * c2.y * w;
        sumw_x+=w;
        sumw_y+=w;
      }
      c.avg_layer += c.layer * w;
      c.z += c2.z * w;
      c.zz += c2.z * c2.z * w;
      sumw += w;
      for (const auto &hit : c2.hits) {
        c.hits.push_back(hit);
        if (hit.layer % 2 == 0) {
          c.strips_evenlayer.push_back(hit.strip);
        } else {
          c.strips_oddlayer.push_back(hit.strip);
        }
        c.strips.push_back(hit.strip);
        s_perlayer.at(hit.layer).push_back(hit.strip);
      }
      c.strips_per_layer = s_perlayer;
    }

    // sort cluster strips                                                                                                                                                                                                                                                                                                                                           
    sort(c.strips_evenlayer.begin(), c.strips_evenlayer.end());
    c.strips_evenlayer.erase(
        unique(c.strips_evenlayer.begin(), c.strips_evenlayer.end()),
        c.strips_evenlayer.end());
    sort(c.strips_oddlayer.begin(), c.strips_oddlayer.end());
    c.strips_oddlayer.erase(
        unique(c.strips_oddlayer.begin(), c.strips_oddlayer.end()),
        c.strips_oddlayer.end());

    // position                                                                                                                                                                                                                                                                                                                                                      
    if(use_toa_){
      c.x /= sumw_x;
      c.y /= sumw_y;
    }
    else{
      c.x /= sumw;
      c.y /= sumw;
    }
    c.avg_layer /= sumw;
    c.z /= sumw;
    c.xx /= sumw;  // now is <x^2>                                                                                                                                                                                                                                                                                                                                   
    c.yy /= sumw;
    c.zz /= sumw;
    c.xx = sqrt(c.xx - c.x * c.x);  // now is sqrt(<x^2>-<x>^2)                                                                                                                                                                                                                                                                                                      
    c.yy = sqrt(c.yy - c.y * c.y);
    c.zz = sqrt(c.zz - c.z * c.z);
  }
  if(debug) {
    cout << "done with 3d rebuilding, size " << clusters3d.size() << endl;
    for (auto &c : clusters3d) c.Print3d();
  }
}
  
std::vector<Cluster> ClusterBuilder::Build2DClustersPerLayer(
    std::vector<Hit> hits) {
  // map hits by id
  std::map<int, Hit> hits_by_id;
  std::map<int, Hit> hits_by_id_temp;
  for (auto &hit : hits) {
    hits_by_id[hit.rawid] = hit;
    hits_by_id_temp[hit.rawid] = hit;
  }

  // 2D algo
  // while (nhits_to_cluster > 0) {
  //  - find the highest energy hit as seed
  //  - loop for neighbor < num_neighbors
  //    - find neighboring strips to the hits in the current cluster
  //      and use those that are below max_xy_2d distance, if using TOA
  //    - erase the hits used from hits_to_cluster
  //  - try to find a next seed if not all the hits have been used
  // }
  
  std::vector<Cluster> clusters;

  int iclus = 0;
  while(!hits_by_id_temp.empty()) {

    // convert map to vector
    std::vector<Hit> dhits;
    MapToVec(hits_by_id_temp, dhits);
      
    // get hit with highest energy as seed
    ESort(dhits);
    auto & seed = dhits.at(0);

    Cluster c;
    if (seed.e > seed_threshold_2d_ ) {
      hits_by_id[seed.rawid].used = true;
      hits_by_id_temp.erase(seed.rawid);
      c.hits.push_back(seed);
      c.e = seed.e;
      c.x = seed.x;
      c.y = seed.y;
      c.z = seed.z;
      c.seed = seed.rawid;
      c.layer = seed.layer;
      c.section = seed.section;
      c.strips.push_back(seed.strip);
    }
    else {
      break;
    }
      
    // add hits to cluster neighbors
    int i_neighbor = 0;    
    while (i_neighbor < num_neighbors_) {
      int num_added = 0;
      for (const auto &hit : c.hits) {
	// if(debug) {
	//   cout << "Looking at neigbors of hit, for the " << i_neighbor<< " time " << endl;
	//   hits_by_id[hit.rawid].Print();
	// }
	for (auto n : geom->strip_neighbors[hit.rawid]) {
	  if (hits_by_id.count(n) &&
	      !hits_by_id[n].used &&
	      hits_by_id[n].e > neighbor_threshold_2d_) {
	    
	    // if using TOA information, use the xy distance
	    // to further prune the neighbors to the seed 2d clusters
	    if (use_toa_ && !isStripNeighbor(hits_by_id[n], hit, max_xy_2d_)) {
	      if (debug) {
		double d = sqrt(pow(hits_by_id[n].x - hit.x, 2) +
                                pow(hits_by_id[n].y - hit.y, 2));
		cout << "  The following hit is discarded because the TOA xy distance is larger than max_xy_2d_: " << max_xy_2d_ << " d " << d << endl;
		hits_by_id[n].Print();
	      }
	      continue;
	    }
	    // else {
	    //   if (debug) {
	    // 	double d = sqrt(pow(hits_by_id[n].x - hit.x, 2) +
	    // 			pow(hits_by_id[n].y - hit.y, 2));
	    // 	cout << "  The following hit is considered because the TOA xy distance is " << d << endl;
	    // 	hits_by_id[n].Print();
	    //   }
	    // }
	    
	    hits_by_id[n].used = true;
	    c.hits.push_back(hits_by_id[n]);

	    hits_by_id_temp.erase(n);
	    num_added++;
	  }
	}
	// if(debug) 
	//   cout << "end hit\n " << endl;
      }
      i_neighbor++;
      if(num_added==0) break;
    }

    iclus++;
    clusters.push_back(c);
  }
  
  // rebuild cluster properties based on the new hits
  // do energy weighting here
  ReBuild2DCluster(clusters);
  if (debug) {
    cout << "\n Preliminary hits and 2d clusters " << endl;
    for (auto &hitpair : hits_by_id) hitpair.second.Print();
    for (auto &c : clusters) c.Print();
  }

  if(is_merge_2d_) {
    if(clusters.size() > 0) {
      // see if its possible to combine clusters
      std::map<int, int> cluster_match;
      for (size_t i = 0; i < clusters.size(); ++i) {
	cluster_match[i] = i;
      }
      
      for (size_t i = 0; i < clusters.size()-1; ++i) {
	for (size_t j = i+1; j < clusters.size(); ++j) {
	  if(i!=j) {
	    double distance = sqrt(pow(clusters.at(i).y - clusters.at(j).y, 2) + pow(clusters.at(i).y - clusters.at(j).y,2) );
	    if(distance <= max_xy_merge_) {
	      cluster_match[j] = cluster_match[i];
	      if(debug) {
		cout << "distance between cluster indices "<< i << " " << j << " is: " << distance << " and below " << max_xy_merge_ << " merging.. " << cluster_match[j] << endl;
	      }
	    }
	    // else {
	    //   cout << "unmerged " << i << " " << cluster_match[i] << " and j " << j << " " << cluster_match[j] << " distance " << distance << endl;
	    // }
	  }
	}
      }
      
      // build a map w cluster index, and hits
      std::vector<int> to_erase;
      for(auto clus : cluster_match) {
	if(clus.second != clus.first) {
	  auto icluster = clusters.at(clus.first);
	  for(auto hit: icluster.hits) 
	    clusters.at(clus.second).hits.push_back(hit);
	  to_erase.push_back(clus.first);
	}
      }
      for(auto index: to_erase)
	clusters.erase(clusters.begin() + index);
      
      // rebuild clusters
      ReBuild2DCluster(clusters);
      if (debug) {
	cout << "\nRe-defined hits and 2d clusters after merging" << endl;
	for (auto &hitpair : hits_by_id) hitpair.second.Print();
	for (auto &c : clusters) c.Print();
      }
    }
  }

  return clusters;
}

void ClusterBuilder::Build2DClusters() {
  // loop over sections
  // TODO: for now only use back hcal
  for (int section = 0; section < 1; section++) {
    // map hits by layer
    std::map<int, std::vector<Hit> > layer_hits;
    for (const auto hit : all_hits) {
      auto l = hit.layer;
      if (hit.section != section) continue;

      // if not using TOA, only consider hits from a given layer parity..
      if(!use_toa_ && (l%2)!=layer_parity_) continue;
      
      if (layer_hits.count(l)) {
        layer_hits[l].push_back(hit);
      } else {
        layer_hits[l] = {hit};
      }
    }

    // run clustering in each layer
    for (auto &pair : layer_hits) {
      if (debug) {
        cout << "----- " << endl;
        cout << "Found " << pair.second.size() << " hits in layer "
             << pair.first << endl;
      }
      auto clus = Build2DClustersPerLayer(pair.second);
      all_clusters.insert(all_clusters.end(), clus.begin(), clus.end());
    }
  }
  
  // sort by energy
  clusters_2d = all_clusters;
  ESort(clusters_2d);
}

void ClusterBuilder::Build3DClusters() {
  if (debug) {
    cout << "\n\n3D: Building 3d clusters" << endl;
  }

  // sort 2d clusters by layer
  std::vector<std::vector<Cluster> > layer_clusters;
  layer_clusters.resize(
      layer_max_);  // NOTE: consider only first 100 2d clusters
  for (auto &clus : all_clusters) {
    layer_clusters[clus.layer].push_back(clus);
  }

  // sort clusters in each layer by energy
  for (auto &clusters : layer_clusters) ESort(clusters);

  if (debug) {
    cout << "\nSorted 2d inputs (in each layer by energy)" << endl;
    for (auto &clusters : layer_clusters)
      for (auto &c : clusters) c.Print();
  }

  // identify cluster seed layer (max shower)
  // repeat, if not all clusters are used
  bool building = true;
  std::vector<Cluster> clusters3d;
  while (building) {
    Cluster cluster3d;
    cluster3d.is2D = false;

    // find seed, from global shower max
    // derived from the first available cluster in each layer
    int layer_showermax = -1;
    double e_showermax = seed_threshold_3d_;
    cout << "finding shower max " << endl;
    for (auto &clusters : layer_clusters) {
      // since the layer clusters are sorted by energy we just check the first element
      if (clusters.size() > 0) {
        if (clusters.at(0).e > e_showermax) {
          layer_showermax = clusters.at(0).layer;
          e_showermax = clusters.at(0).e;
        }
      }
    }

    if (layer_showermax == -1) break;

    auto &clusters2d = layer_clusters[layer_showermax];
    if (debug) {
      cout << "\nFound showermax: ";
      clusters2d.at(0).Print();
    }

    cluster3d.clusters2d = {clusters2d.at(0)};
    cluster3d.depth = 1;
    cluster3d.first_layer = layer_showermax;
    cluster3d.last_layer = layer_showermax;
    cluster3d.avg_layer = layer_showermax;
    cluster3d.seed = clusters2d.at(0).seed;
    clusters2d.erase(clusters2d.begin());
    
    // loop to front of hcal from the showermax, and then to back of hcal
    int found_in_last_layer = 1;
    for (int ilayer = 1; ilayer < layer_max_; ilayer++) {
      int test_layer = ilayer;
      if (ilayer < layer_max_ - layer_showermax) {
        test_layer = layer_showermax + ilayer;
      } else {
        test_layer = layer_max_ - ilayer - 1;
      }

      // get all possible 2d clusters from that layer
      // and only look into layers that have 2d clusters
      auto &clusters2d = layer_clusters[test_layer];
      if (clusters2d.size() <= 0) {
        found_in_last_layer = 0;
        continue;
      }
      
      // if layer is shower max make sure you do not include the seed
      if(test_layer == layer_showermax) {
	if (debug) {
	  cout << "Same layer as shower max " << clusters2d.size() << endl;
	}
      }

      // get last 2d seed from 3d cluster
      auto &last_cluster = cluster3d.clusters2d.back();
      int last_cluster_layer = last_cluster.layer;
      auto &last_seed2d = last_cluster.seed;
      
      // TODO:
      // sometimes layers next to the shower max might not have hits
      // so we modify the last seed reference here
      if ((test_layer <= layer_showermax - 1) &&
          (last_cluster.layer > layer_showermax) && found_in_last_layer == 0) {
	last_cluster_layer = cluster3d.clusters2d.front().layer;
	last_seed2d = cluster3d.clusters2d.front().seed;
        if (debug) {
          cout << "did not find in last layer ";
          cout << " new seed layer is " << last_cluster_layer << endl;
        }
      }

      if (debug) {
        cout << "Look at " << clusters2d.size() << " clusters in layer "
             << test_layer;
        cout << " with last seed layer " << last_cluster_layer << endl;
      }

	
      int found_this_layer = 0;
      for (int iclus2d = 0; iclus2d < clusters2d.size(); iclus2d++) {
        // require an energy threshold for the neighboring 2d clusters
        auto &seed2d = clusters2d[iclus2d].seed;

        if (clusters2d[iclus2d].e <= neighbor_threshold_3d_) {
          if (debug) {
            cout << " Cluster does not pass neighbor threshold " << endl;
          }
          continue;
        }

        // sometimes layers next to the layer of the last seed might not have
        // hits we can have a threshold of 1 (maybe? make it configurable)
        if (found_in_last_layer == 0 &&
            ((last_cluster_layer < test_layer - close_layer_) ||
             (last_cluster_layer > test_layer + close_layer_))) {
          if (debug) {
            cout << " Cluster is not a close enough neighbor " << last_cluster_layer << " "
                 << test_layer << endl;
          }
          continue;
        }
        if (found_in_last_layer == 1 &&
            !geom->CheckLayerNeighbor(last_seed2d, seed2d) &&
            last_seed2d != seed2d) {
          if (debug) {
            cout << " Cluster is not a layer neighbor " << last_cluster_layer
                 << " " << test_layer << endl;
          }
          continue;
        }

        if (use_toa_) {
          // check distance between cluster centers
          Hit clus2d_center;
          clus2d_center.x = clusters2d[iclus2d].x;
          clus2d_center.y = clusters2d[iclus2d].y;
          Hit lastclus2d_center;
          lastclus2d_center.x = last_cluster.x;
          lastclus2d_center.y = last_cluster.y;
          if (!isStripNeighbor(clus2d_center, lastclus2d_center, max_xy_3d_)) {
            if (debug) {
              cout << " The following 2D cluster is not a neighbor: ";
              clusters2d[iclus2d].Print();
            }
            continue;
          }
        }

        if (debug) {
          cout << " Extend with this cluster: ";
          clusters2d[iclus2d].Print();
        }

        if (test_layer < cluster3d.first_layer)
          cluster3d.first_layer = test_layer;
        if (test_layer > cluster3d.last_layer)
          cluster3d.last_layer = test_layer;
	
        cluster3d.clusters2d.push_back(clusters2d[iclus2d]);
        cluster3d.depth++;
        clusters2d.erase(clusters2d.begin() + iclus2d);

        // once it found a neighboring 2d cluster, proceed to next layer
        // TODO: are all the 3D clusters associated to only 1 2d cluster in each
        // layer?? perhaps not..
        found_this_layer = 1;
        break;
      }
      if (found_this_layer == 0) found_in_last_layer = 0;

    }  // end loop over layers

    if(debug) {
      cout << "done with 3d cluster and its depth is " << cluster3d.depth << endl;
      cout << " The following are 2d clusters that are left: " << endl;
      for (auto &clusters : layer_clusters)
	for (auto &c : clusters) c.Print();
      cout << "\n" << endl;
    }
    if (cluster3d.depth == 0) {
      cout << "building false " << endl;
      building = false;
    } else {
      if (cluster3d.depth >= 1) clusters3d.push_back(cluster3d);
    }
  }  // building

  if(debug)
    cout << " done with building " << endl;
  
  // post-process
  ReBuild3DCluster(clusters3d);

  if (debug) {
    cout << "\nAfter postprocess, found 3d clusters" << endl;
    for (auto &c : clusters3d) c.Print3d();
  }

  if(clusters3d.size() > 0) {
    Merge3DClusters(clusters3d);
  }
  
  all_clusters.clear();
  all_clusters.insert(all_clusters.begin(), clusters3d.begin(),
                      clusters3d.end());
}

void ClusterBuilder::Merge3DClusters(std::vector<Cluster> & clusters3d) {
  cout << "size " <<  clusters3d.size() << endl;
  std::map<int, int> cluster_match;
  for (size_t i = 0; i < clusters3d.size(); ++i) {
     cluster_match[i] = i;
  }

  for (size_t i = 0; i < clusters3d.size()-1; ++i) {
    for (size_t j = i+1; j < clusters3d.size(); ++j) {
      if(i!=j) {
	// x and y distance
	double distance_xy = sqrt(pow(clusters3d.at(i).x - clusters3d.at(j).x, 2) + pow(clusters3d.at(i).y - clusters3d.at(j).y,2) );
	//cout << "d " << distance_xy << endl;
	// z distance
	if(distance_xy <= max_xy_merge_) {
	  cluster_match[j] = cluster_match[i];
	  cout << "distance between cluster indices "<< i << " " << j << " is: " << distance_xy << " and below " << max_xy_merge_ << " merging.. " << cluster_match[j] << endl;
	}
	else{
	  cout << "unmerged " << i << " " << cluster_match[i] << " and j " << j << " " << cluster_match[j] << " distance " << distance_xy << endl;
	}
      }
    }
  }

  std::vector<int> to_erase;
  for(auto clus : cluster_match) {
    cout << "cluster match " << clus.first << " " << clus.second << endl;
    if(clus.second != clus.first) {
      auto icluster = clusters3d.at(clus.first);
      for(auto cluster2d: icluster.clusters2d)
	clusters3d.at(clus.second).clusters2d.push_back(cluster2d);
      //for(auto hit: icluster.hits)
      //	clusters3d.at(clus.second).hits.push_back(hit);
      cout << "erasing index " << clus.first << endl;
      to_erase.push_back(clus.first);
    }
  }                                                                                                                                                                                                                                                                           
  for(auto index: to_erase)
    clusters3d.erase(clusters3d.begin() + index);

  ReBuild3DCluster(clusters3d);

  if (debug) {
    cout << "\nRe-defined hits and 3d clusters " << endl;
    for (auto &c : clusters3d) c.Print3d();
  }                                    
  
}

void ClusterBuilder::BuildClusters() {
  if (debug) {
    cout << "--------\nAll hits" << endl;
    for (auto &hit : all_hits) hit.Print();
  }

  // cluster the hits in each layer
  Build2DClusters();

  // cluster the hits across layers
  Build3DClusters();

  // sort by energy
  ESort(all_clusters);
}

}  // namespace hcal