/**
 * @file HcalVetoResult.h
 * @brief Class used to encapsulate the results obtained from
 *        HcalVetoProcessor.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef HCAL_EVENT_HCALVETORESULT_H_
#define HCAL_EVENT_HCALVETORESULT_H_

#include <fire/io/Access.h>

//----------//
//   ROOT   //
//----------//
#include "TObject.h"  //For ClassDef

//----------//
//   LDMX   //
//----------//
#include "Hcal/Event/HcalHit.h"

namespace ldmx {

class HcalVetoResult {
 public:
  /** Constructor */
  HcalVetoResult();

  /** Destructor */
  ~HcalVetoResult();

  /** Reset the object. */
  void clear();

  /** Print out the object */
  void Print() const;

  /** Checks if the event passes the Hcal veto. */
  bool passesVeto() const { return passesVeto_; };

  /** @return The maximum PE HcalHit. */
  inline ldmx::HcalHit getMaxPEHit() const { return maxPEHit_; }

  /**
   * Sets whether the Hcal veto was passed or not.
   *
   * @param passesVeto Veto result.
   */
  inline void setVetoResult(const bool& passesVeto = true) {
    passesVeto_ = passesVeto;
  }

  /**
   * Set the maximum PE hit.
   *
   * @param maxPEHit The maximum PE HcalHit
   */
  inline void setMaxPEHit(const ldmx::HcalHit maxPEHit) {
    maxPEHit_ = maxPEHit;
  }

 private:
  /** Reference to max PE hit. */
  ldmx::HcalHit maxPEHit_;

  /** Flag indicating whether the event passes the Hcal veto. */
  bool passesVeto_{false};

  friend class fire::io::access;
  template<typename Data>
  void attach(Data& d) {
    d.attach("maxPEHit", maxPEHit_);
    d.attach("passesVeto", passesVeto_);
  }

  ClassDef(HcalVetoResult, 2);

};  // HcalVetoResult
}  // namespace ldmx

#endif  // HCAL_EVENT_HCALVETORESULT_H_
