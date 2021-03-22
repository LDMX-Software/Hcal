#include "Hcal/Event/HcalMIPTracks.h"

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>

//-------------//
//   ldmx-sw   //
//-------------//
#include "Hcal/Event/HcalHit.h"

ClassImp(ldmx::HcalMIPTracks)

    namespace ldmx {
  HcalMIPTracks::HcalMIPTracks() {}

  HcalMIPTracks::~HcalMIPTracks() {}

  void HcalMIPTracks::Clear() {
  }

  void HcalMIPTracks::Print() const {
    std::cout << "[ HcalMIPTracks ]: Track : "
              << std::endl;
  }
}  // namespace ldmx
