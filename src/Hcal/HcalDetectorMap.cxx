#include "Hcal/HcalDetectorMap.h"

#include <sstream>

#include "Framework/ConditionsObjectProvider.h"
#include "Framework/EventHeader.h"

namespace hcal {

class HcalDetectorMapLoader : public framework::ConditionsObjectProvider {
 public:
  HcalDetectorMapLoader(const std::string& name, const std::string& tagname,
                        const framework::config::Parameters& parameters,
                        framework::Process& process)
      : ConditionsObjectProvider(HcalDetectorMap::CONDITIONS_OBJECT_NAME,
                                 tagname, parameters, process),
        the_map_{nullptr} {
    want_d2e_ = parameters.getParameter<bool>("want_d2e");
    connections_table_ = parameters.getParameter<std::string>("connections_table");
  }

  virtual std::pair<const framework::ConditionsObject*,
                    framework::ConditionsIOV>
  getCondition(const ldmx::EventHeader& context) {
    if (!the_map_) {
      the_map_ = new HcalDetectorMap(connections_table_, want_d2e_);
    }

    return std::make_pair(
        the_map_, framework::ConditionsIOV(context.getRun(), context.getRun(),
                                          true, true));
  }

  /**
   * Take no action on release, as the object is permanently owned by the
   * Provider
   */
  virtual void releaseConditionsObject(const framework::ConditionsObject* co) {}

 private:
  HcalDetectorMap* the_map_;
  std::string connections_table_;
  bool want_d2e_;
};

HcalDetectorMap::HcalDetectorMap(const std::string& connections_table, bool want_d2e)
    : framework::ConditionsObject(CONDITIONS_OBJECT_NAME),
      ldmx::ElectronicsMap<ldmx::HcalElectronicsID, ldmx::HcalDigiID>(want_d2e) {
  
  this->clear();
  conditions::StreamCSVLoader csv(connections_table);
  while (csv_loader.nextRow()) {
    /** Column Names
     * "HGCROC" "Channel" "CMB" "Quadbar" "Bar" "Plane"
     *
     * Not confident about this!!!
     */
    ldmx::HcalDigiID detid(0 /*section - only one section during test beam*/, 
        csv.getInteger("Plane") /*layer*/,
        csv.getInteger("Bar") /*strip*/,
        csv.getInteger("Quadbar")-1 /*end??*/); //Quadbar given as 1 or 2
    ldmx::HcalElectronicsID eleid(
        0 /*fpga - only one FPGA during test beam*/
        csv.getInteger("HGCROC") /*fiber*/,
        csv.getInteger("Channel") /*channel*/);

     if (this->exists(eleid)) {
       std::stringstream ss;
       ss << "Two different mappings for electronics channel " << eleid;
       EXCEPTION_RAISE("DuplicateMapping", ss.str());
     }
     this->addEntry(eleid, detid);
  }
}

}  // namespace hcal
DECLARE_CONDITIONS_PROVIDER_NS(hcal, HcalDetectorMapLoader);
