from LDMX.Framework import ldmxcfg

class HcalNewClusterProducer(ldmxcfg.Producer) :
    """Configuration for cluster producer in the HCal

    Examples
    --------
        import LDMX.Hcal.cluster as hcal_cluster
        p.sequence.append( hcal_cluster.HcalNewClusterProducer() ) 
    """

    def __init__(self, instance_name = 'hcalClusters',
                 pass_name = '', coll_name = 'HcalRecHits',
                 cluster2d_coll_name = 'Hcal2DClusters',
                 cluster3d_coll_name = 'Hcal3DClusters'):
        super().__init__(instance_name,'hcal::HcalNewClusterProducer','Hcal')

        self.coll_name = coll_name
        self.pass_name = pass_name
        self.cluster2d_coll_name = cluster2d_coll_name
        self.cluster3d_coll_name = cluster3d_coll_name

        # energy thresholds
        self.noise_threshold = 0.01 # MeV
        self.seed_threshold_2d = 0.1
        self.neighbor_threshold_2d = 0.01

        self.seed_threshold_3d = 4 # MeV ( a MIP?)
        self.neighbor_threshold_2d = 1
        
        # number of neighboring strips
        self.num_neighbors = 4

        # max xy (mm)
        self.max_xy_2d = 3*50.
        self.max_xy_3d = 8*50.

        self.use_toa = True
