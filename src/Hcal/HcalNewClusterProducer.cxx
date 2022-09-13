#include "Hcal/HcalNewClusterProducer.h"

#include "Hcal/ClusterBuilder.h"
#include "Hcal/HcalReconConditions.h"

namespace hcal {

void HcalNewClusterProducer::configure(framework::config::Parameters& p) {
  pass_name_ = p.getParameter("pass_name", pass_name_);
  coll_name_ = p.getParameter("coll_name", coll_name_);

  cluster2d_coll_name_ =
      p.getParameter("cluster2d_coll_name", cluster2d_coll_name_);
  cluster3d_coll_name_ =
      p.getParameter("cluster3d_coll_name", cluster3d_coll_name_);

  num_neighbors_ = p.getParameter<int>("num_neighbors");
  noise_threshold_ = p.getParameter<double>("noise_threshold");
  seed_threshold_2d_ = p.getParameter<double>("seed_threshold_2d");
  neighbor_threshold_2d_ = p.getParameter<double>("neighbor_threshold_2d");
  seed_threshold_3d_ = p.getParameter<double>("seed_threshold_3d");
  neighbor_threshold_3d_ = p.getParameter<double>("neighbor_threshold_3d");
  max_xy_2d_ = p.getParameter<double>("max_xy_2d");
  max_xy_3d_ = p.getParameter<double>("max_xy_3d");
  use_toa_ = p.getParameter<bool>("use_toa");
}

void HcalNewClusterProducer::produce(framework::Event& event) {
  const auto& hcalGeometry = getCondition<ldmx::HcalGeometry>(
      ldmx::HcalGeometry::CONDITIONS_OBJECT_NAME);

  auto hcalRecHits = event.getCollection<ldmx::HcalHit>(coll_name_, pass_name_);

  ClusterGeometry clusterGeometry(hcalGeometry);

  ClusterBuilder builder;
  builder.SetThresholds2D(seed_threshold_2d_, neighbor_threshold_2d_);
  builder.SetThresholds3D(seed_threshold_3d_, neighbor_threshold_3d_);
  builder.SetNeighbors(num_neighbors_);
  builder.SetMaxXY(max_xy_2d_, max_xy_3d_);
  builder.SetTOA(use_toa_);
  builder.SetClusterGeo(&clusterGeometry);
  for (auto const& h : hcalRecHits) {
    // quality cuts for hits entering clustering
    if (h.getEnergy() < noise_threshold_) continue;
    builder.AddHit(h);
  }
  builder.BuildClusters();

  auto clusters_2d = builder.Get2DClusters();
  auto clusters_3d = builder.Get3DClusters();

  std::vector<ldmx::HcalCluster> hcalClusters_2d, hcalClusters_3d;
  for (const auto& c : clusters_2d) {
    ldmx::HcalCluster cluster;
    cluster.setEnergy(c.e);
    cluster.setNHits(c.hits.size());
    cluster.setCentroidXYZ(c.x, c.y, c.z);
    cluster.setRMSXYZ(c.xx, c.yy, c.zz);
    cluster.addStrips(c.strips);

    cluster.setSeedEnergy(c.hits.at(0).e);
    cluster.setLayer(c.layer);

    hcalClusters_2d.push_back(cluster);
  }

  // add collection to event bus
  event.add(cluster2d_coll_name_, hcalClusters_2d);

  for (const auto& c : clusters_3d) {
    ldmx::HcalCluster cluster;
    cluster.setEnergy(c.e);
    cluster.setNHits(c.hits.size());
    cluster.setCentroidXYZ(c.x, c.y, c.z);
    cluster.setRMSXYZ(c.xx, c.yy, c.zz);
    cluster.addStrips(c.strips);

    cluster.setStripsOdd(c.strips_oddlayer);
    cluster.setStripsEven(c.strips_evenlayer);

    cluster.setDepth(c.depth);
    cluster.setSeedEnergy(c.clusters2d.at(0).e);
    cluster.setN2DClusters(c.clusters2d.size());
    cluster.setLayer(c.first_layer);

    hcalClusters_3d.push_back(cluster);
  }
  // add collection to event bus
  event.add(cluster3d_coll_name_, hcalClusters_3d);
}

}  // namespace hcal
DECLARE_PRODUCER_NS(hcal, HcalNewClusterProducer);
