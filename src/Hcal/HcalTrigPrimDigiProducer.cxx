#include "Hcal/HcalTrigPrimDigiProducer.h"

#include "DetDescr/HcalGeometry.h"
#include "Hcal/HcalTriggerGeometry.h"
#include "Recon/Event/CaloTrigPrim.h"
#include "Recon/Event/HgcrocDigiCollection.h"
#include "Recon/Event/HgcrocTrigDigi.h"
#include "Tools/HgcrocTriggerCalculations.h"

namespace hcal {
HcalTrigPrimDigiProducer::HcalTrigPrimDigiProducer(const std::string& name,
                                                   framework::Process& process)
    : Producer(name, process) {}

void HcalTrigPrimDigiProducer::configure(framework::config::Parameters& ps) {
  digiCollName_ = ps.getParameter<std::string>("digiCollName");
  digiPassName_ = ps.getParameter<std::string>("digiPassName");
  condObjName_ =
      ps.getParameter<std::string>("condObjName", "HcalTrigPrimDigiConditions");
}

void HcalTrigPrimDigiProducer::produce(framework::Event& event) {
  const HcalTriggerGeometry& geom = getCondition<HcalTriggerGeometry>(
      HcalTriggerGeometry::CONDITIONS_OBJECT_NAME);

  const ldmx::HgcrocDigiCollection& hcalDigis =
      event.getObject<ldmx::HgcrocDigiCollection>(digiCollName_, digiPassName_);

  // get the calibration object
  const conditions::IntegerTableCondition& conditions =
      getCondition<conditions::IntegerTableCondition>(condObjName_);

  // construct the calculator...
  ldmx::HgcrocTriggerCalculations calc(conditions);

  bool debug=false;
  if(debug) std::cout << "\nNewEvent\n";

  // TMP for cross-check only
  ldmx::CaloTrigPrimCollection debugDaqDigis;
  
  // Loop over the digis
  for (unsigned int ix = 0; ix < hcalDigis.getNumDigis(); ix++) {
    const ldmx::HgcrocDigiCollection::HgcrocDigi pdigi = hcalDigis.getDigi(ix);
    ldmx::HcalTriggerID tid = geom.belongsToQuad(ldmx::HcalDigiID(pdigi.id()));

    if (!tid.null()) {
      int tot = 0;
      if (pdigi.soi().isTOTComplete()) tot = pdigi.soi().tot();
      calc.addDigi(pdigi.id(), tid.raw(), pdigi.soi().adc_t(), tot);
      debugDaqDigis.push_back( ldmx::CaloTrigPrim(pdigi.id(), pdigi.soi().adc_t()));
      if(debug) std::cout << "DAQ digi with ADC = " << pdigi.soi().adc_t() << std::endl;
    }
  }
  event.add("debugHcalDaqDigis", debugDaqDigis);

  // Now, we compress the digis
  calc.compressDigis(4);
  const float hgc_compress_factor = 2;

  // const std::map<unsigned int, uint8_t> results;
  const std::map<unsigned int, uint8_t>& results = calc.compressedEnergies();
  ldmx::HgcrocTrigDigiCollection tdigis;
  ldmx::CaloTrigPrimCollection tdigisUC; // sums w/o compression applied

  for (auto result : results) {
    if (result.second > 0) {
      int uncompress = hgc_compress_factor*ldmx::HgcrocTrigDigi::compressed2Linear(result.second);
      if(debug) std::cout << "TRIG digi with ADC = " << uncompress
                          << " (compressed: "<<(unsigned int)(result.second)<<")" << std::endl;
      tdigis.push_back(ldmx::HgcrocTrigDigi(result.first, result.second));
      tdigisUC.push_back( ldmx::CaloTrigPrim(result.first, uncompress));
    }
  }

  // Form TriggerQuads from 4 CMS quads
  //   (TODO: and compressed their energies)
  trigQuadMap.clear();
  for (auto result : results) {
    if (result.second > 0) {
      const ldmx::HcalTriggerID quad_id(result.first);
      const std::vector<ldmx::HcalDigiID> precisions_ids =
          geom.contentsOfQuad(quad_id);
      if (precisions_ids.size() == 0) {
        std::cout << "Failing HcalTriggerID(" << quad_id.section() << ','
                  << quad_id.layer() << ',' << quad_id.superstrip() << ','
                  << quad_id.end() << ')' << std::endl;
        EXCEPTION_RAISE("TrigToPrecIDMismatch",
                        "Attempted to lookup a nonexistent HcalDigiID from an "
                        "HcalTriggerID");
      }
      const ldmx::HcalDigiID prec_id(precisions_ids.front().raw());
      const ldmx::HcalTriggerID trigQuad_id = geom.belongsToTrigQuad(prec_id);
      auto ptr = trigQuadMap.find(trigQuad_id.raw());
      int linear_charge =
          hgc_compress_factor *
          ldmx::HgcrocTrigDigi::compressed2Linear(result.second);
      if (ptr != trigQuadMap.end()) {
        ptr->second += linear_charge;
      } else {
        trigQuadMap[trigQuad_id.raw()] = linear_charge;
      }
    }
  }

  ldmx::CaloTrigPrimCollection trigQuads;
  for (auto result : trigQuadMap) {
    if (result.second > 0) {
      if(debug) std::cout << "TRIG quad with ADC = " << result.second << std::endl;
      trigQuads.push_back(ldmx::CaloTrigPrim(result.first, result.second));
    }
  }

  event.add(getName(), tdigis); // hcalTrigPrimDigis
  event.add(getName()+"Uncompress", tdigisUC);
  event.add("hcalOneEndedTrigQuads", trigQuads);
}
}  // namespace hcal
DECLARE_PRODUCER_NS(hcal, HcalTrigPrimDigiProducer);
