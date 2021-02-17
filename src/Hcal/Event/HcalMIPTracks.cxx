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
      //trackHits_;
      x_ = -9999.;
      y_ = -9999.;
      dx_ = -9999.;
      dy_ = -9999.;
      xx_ = -9999.;
      xy_ = -9999.;
      xdx_ = -9999.;
      xdy_ = -9999.;
      yy_ = -9999.;
      ydx_ = -9999.;
      ydy_ = -9999.;
      dxdx_ = -9999.;
      dxdy_ = -9999.;
      dydy_ = -9999.;
      nTracks_ = 0;
      isTriggered_ = false;
  }

  void HcalMIPTracks::Print() const {
    std::cout << "[ HcalMIPTracks ]: Track : "
              << std::endl;
  }
}  // namespace ldmx
