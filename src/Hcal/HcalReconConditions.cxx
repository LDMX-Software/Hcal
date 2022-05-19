#include "Hcal/HcalReconConditions.h"

#include "Framework/ConditionsObjectProvider.h"
#include "Framework/EventHeader.h"

namespace hcal {

const std::string HcalReconConditions::CONDITIONS_NAME = "HcalReconConditions";

HcalReconConditions::HcalReconConditions(const conditions::DoubleTableCondition& adc_ped, 
    const conditions::DoubleTableCondition& adc_gain,
    const conditions::DoubleTableCondition& tot_ped,
    const conditions::DoubleTableCondition& tot_gain)
  : framework::ConditionsObject(HcalReconConditions::CONDITIONS_NAME),
    adc_pedestals_{adc_ped}, adc_gains_{adc_gain},
    tot_pedestals_{tot_ped}, tot_gains_{tot_gain} {}

/**
 * a helpful interface for grabbing the parent conditions at once
 *
 * the most complicated task it does is calculate the "minimum" IOV from 
 * the parent conditions to make sure that it is updated as soon as it needs
 * to be updated.
 *
 * @TODO Right now, no minimum IOV calculation is being done. We just have
 * our IOV be one run at a time.
 */
class HcalReconConditionsProvider : public framework::ConditionsObjectProvider {
  /// name of condition object for hcal adc gains
  std::string adc_gain_;
  /// name of condition object for hcal adc pedestals
  std::string adc_ped_;
  /// name of condition object for hcal tot gains
  std::string tot_gain_;
  /// name of condition object for hcal tot pedestals
  std::string tot_ped_;
 public:
  /**
   * Retrieve the name of the parent conditions from the configuration
   *
   * @throw Exception if name is not what it should be
   *
   * @param[in] name  should be HcalReconConditions::CONDITIONS_NAME
   * @param[in] tagname
   * @param[in] parameters python configuration parameters
   * @param[in] proc handle to current process
   */
  HcalReconConditionsProvider(const std::string& name, const std::string& tagname,
                              const framework::config::Parameters& parameters,
                              framework::Process& proc)
    : ConditionsObjectProvider(hcal::HcalReconConditions::CONDITIONS_NAME,
                               tagname, parameters, proc) {
      if (name != HcalReconConditions::CONDITIONS_NAME) {
        EXCEPTION_RAISE("BadConfig",
            "The name provided to HcalReconConditionsProvider "+name
            +" is not equal to the expected name "+HcalReconConditions::CONDITIONS_NAME);
      }
      adc_gain_ = parameters.getParameter<std::string>("adc_gain");
      adc_ped_ = parameters.getParameter<std::string>("adc_ped");
      tot_gain_ = parameters.getParameter<std::string>("tot_gain");
      tot_ped_ = parameters.getParameter<std::string>("tot_ped");
    }

  /**
   * Get the wrapped condition
   *
   * This is where we deduce the "minimum" IOV.
   * @note Right now, we just have the IOV go one run at a time.
   *
   * @note Expects the parent condition tables to all be conditions::DoubleTableCondition
   *
   * @see requestParentCondition for how we get the parent condition tables
   */
  virtual std::pair<const framework::ConditionsObject*, framework::ConditionsIOV>
  getCondition(const ldmx::EventHeader& context) final override {
    // requestParentCondition does check current context for validity 
    // to avoid extra constructions
    auto [ adc_gain_co, adc_gain_iov ] = requestParentCondition(adc_gain_, context);
    auto [ adc_ped_co , adc_ped_iov  ] = requestParentCondition(adc_ped_ , context);
    auto [ tot_gain_co, tot_gain_iov ] = requestParentCondition(tot_gain_, context);
    auto [ tot_ped_co , tot_ped_iov  ] = requestParentCondition(tot_ped_ , context);
    
    // deduce "minimum" IOV
    //  Framework #56 : https://github.com/LDMX-Software/Framework/issues/56
    auto min_iov = adc_ped_iov;
    // use std::move(min_iov) in return statement below
    
    // wrap
    framework::ConditionsObject* co = new hcal::HcalReconConditions(
        dynamic_cast<const conditions::DoubleTableCondition&>(*adc_ped_co),
        dynamic_cast<const conditions::DoubleTableCondition&>(*adc_gain_co),
        dynamic_cast<const conditions::DoubleTableCondition&>(*tot_ped_co),
        dynamic_cast<const conditions::DoubleTableCondition&>(*tot_gain_co)
        );

    return { co , framework::ConditionsIOV(context.getRun(),context.getRun()) };
  }
};

}  // namespace hcal

DECLARE_CONDITIONS_PROVIDER_NS(hcal, HcalReconConditionsProvider);
