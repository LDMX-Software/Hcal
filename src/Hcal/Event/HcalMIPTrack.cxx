#include "Hcal/Event/HcalMIPTrack.h"

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>

//-------------//
//   ldmx-sw   //
//-------------//
#include "Hcal/Event/HcalHit.h"

ClassImp(ldmx::HcalMIPTrack)

    namespace ldmx {
  HcalMIPTrack::HcalMIPTrack() {}

  HcalMIPTrack::~HcalMIPTrack() {}

  void HcalMIPTrack::Clear() {
  }

  void HcalMIPTrack::Print() const {
    std::cout << "[ HcalMIPTrack ]: Track : "
              << std::endl;
  }
}  // namespace ldmx
