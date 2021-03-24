/*
   TemplatedClusterFinder
*/

#ifndef HCAL_TEMPLATEDCLUSTERFINDER_H_ 
#define HCAL_TEMPLATEDCLUSTERFINDER_H_

#include "Hcal/WorkingCluster.h"
#include "TH2F.h"

#include <math.h>
#include <map>

namespace hcal {

    template <class WeightClass>

    class TemplatedClusterFinder {
    
        public:

            void add(const ldmx::HcalHit* eh, const ldmx::HcalGeometry& hex) {
                //This makes a new cluster for every hit... the trick comes in combining the clusters.
                clusters_.push_back(WorkingCluster(eh, hex));
            }
            
    

            static bool compClusters(const WorkingCluster& a, const WorkingCluster& b) {
                return a.centroid().E() > b.centroid().E();
            }

            void cluster(double seed_threshold, double cutoff, double deltaTime) {
                int ncluster = clusters_.size();
                double minwgt = cutoff;
                //std::cout<<"cut off "<<cutoff<<std::endl;
                //sort clusters by energy (at this point 
                std::sort(clusters_.begin(), clusters_.end(), compClusters);
                do {
                    bool any = false;
                    unsigned int mi(0),mj(0);
   
                    int nseeds = 0;
                    // loop over all clusters:
                    for (unsigned int i = 0; i < clusters_.size(); i++) {
                      //                      std::cout<<" i cluster "<<i<<" of "<< clusters_.size()<<" wgt "<<minwgt<< " nclusters " << ncluster<< std::endl;
                        // skip if empty
                        if (clusters_[i].empty()) continue;
                        // check if cluster might be a seed minimum seed:
                        bool iseed = (clusters_[i].centroid().E() >= seed_threshold);
                        if (iseed) {
                            nseeds++;
                            
                        } else {
                            // Since we sorted initially if we find a hit below seed threshold
                            // it is guaranteed that there will be no seeds after.
                            break;
                        }
                        //loop over the rest of the clusters:
                        for (unsigned int j = i + 1; j < clusters_.size(); j++) {
                          //    std::cout<<" j cluster "<<j<<" of "<< clusters_.size()<<std::endl;
                            if (clusters_[j].empty() or (!iseed and clusters_[j].centroid().E() < seed_threshold)) continue;
                            // calculate weights between the two clusters:
                            double wgt = wgt_(clusters_[i],clusters_[j]);// TODO
                            if (!any or wgt < minwgt) { //minwgt begins at cut-off, this probably wants to tell us if these can physically be part of same shower (deltaT, and R cuts)
                                any = true;
                                minwgt = wgt;
                                mi = i;
                                mj = j;
                            }
                            //std::cout<<" end wgt "<<minwgt<<std::endl;
                        }
                    }
                    //if(abs(clusters_[mi].GetTime() - clusters_[mj].GetTime()) > deltaTime) continue;

                    //                    std::cout << clusters_[mi] << " + " << clusters_[mj] << std::endl;
                    nseeds_ = nseeds;
                    transitionWeights_.insert(std::pair<int, double>(ncluster, minwgt));
                    //std::cout<<"Minium Weight "<<minwgt<<std::endl;
                    if (any and minwgt < cutoff) {
                    
                        // put the bigger one in mi
                        if (clusters_[mi].centroid().E() < clusters_[mj].centroid().E()) { std::swap(mi,mj); }
                        // now we have the smallest, merge
                        clusters_[mi].add(clusters_[mj]);
                        clusters_[mj].clear();
                        // decrement cluster count
                        ncluster--;
                    } 
    
                } while (minwgt < cutoff and ncluster > 1);
                finalwgt_ = minwgt;
            }

            double getYMax() const { return finalwgt_; }

            int getNSeeds() const { return nseeds_; }

            std::map<int, double> getWeights() const { return transitionWeights_; }

            std::vector<WorkingCluster> getClusters() const {
                return clusters_;
            }
    
        private:
    
            WeightClass wgt_;
            double finalwgt_;
            int nseeds_;
            std::map<int, double> transitionWeights_;
            std::vector<WorkingCluster> clusters_;
    };
}

#endif
