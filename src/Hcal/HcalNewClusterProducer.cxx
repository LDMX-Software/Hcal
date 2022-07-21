#include "Hcal/HcalNewClusterProducer.h"
#include "Hcal/ClusterBuilder.h"
#include "Hcal/HcalReconConditions.h"

namespace hcal {

void HcalNewClusterProducer::configure(framework::config::Parameters& p) {
  pass_name_ = p.getParameter("pass_name", pass_name_);
  coll_name_ = p.getParameter("coll_name", coll_name_);

  cluster2d_coll_name_ = p.getParameter("cluster2d_coll_name", cluster2d_coll_name_);
  cluster3d_coll_name_ = p.getParameter("cluster3d_coll_name", cluster3d_coll_name_);

  seed_threshold_ = p.getParameter("seed_threshold", seed_threshold_);
  neighbor_threshold_ = p.getParameter("neighbor_threshold", neighbor_threshold_);
  num_neighbors_ = p.getParameter("num_neighbors", num_neighbors_);
  
}

void HcalNewClusterProducer::produce(framework::Event& event) {
  const auto& hcalGeometry = getCondition<ldmx::HcalGeometry>(
      ldmx::HcalGeometry::CONDITIONS_OBJECT_NAME);

  const auto& conditions{
      getCondition<HcalReconConditions>(HcalReconConditions::CONDITIONS_NAME)};

  auto hcalRecHits =
      event.getCollection<ldmx::HcalHit>(coll_name_, pass_name_);

  ClusterGeometry clusterGeometry(hcalGeometry.getStripPositionMap());
  
  ClusterBuilder builder;
  builder.SetThresholds( seed_threshold_, neighbor_threshold_);
  builder.SetNeighbors( num_neighbors_ );
  builder.SetClusterGeo( &clusterGeometry );
  for (auto const& h : hcalRecHits) builder.AddHit(h);
  builder.BuildClusters();
  auto clusters_2d = builder.Get2DClusters();
  //auto clusters_3d = builder.Get3DClusters();
  
  std::vector<ldmx::HcalCluster> hcalClusters_2d, hcalClusters_3d;
  for (const auto &c : clusters_2d) {
    ldmx::HcalCluster cluster;
    cluster.setEnergy(c.e);
    cluster.setNHits(c.hits.size());
    cluster.setCentroidXYZ(c.x,c.y,c.z);
    cluster.setRMSXYZ(c.xx,c.yy,c.zz);
    cluster.addStrips(c.strips);
    cluster.setLayer(c.layer);
    hcalClusters_2d.push_back(cluster);
  }
  
  // add collection to event bus
  event.add(cluster2d_coll_name_, hcalClusters_2d);

  // for (const auto &c : clusters_3d) {
  //   ldmx::HcalCluster cluster;
  //   cluster.setEnergy(c.e);
  //   cluster.setNHits(c.hits.size());
  //   cluster.setCentroidXYZ(c.x,c.y,c.z);
  //   cluster.setRMSXYZ(c.xx,c.yy,c.zz);
  //   cluster.addStrips(c.strips);
  //   hcalClusters_3d.push_back(cluster);
  // }
  // // add collection to event bus
  // event.add(cluster3d_coll_name_, hcalClusters_3d);

}

}  // namespace hcal
DECLARE_PRODUCER_NS(hcal, HcalNewClusterProducer);
