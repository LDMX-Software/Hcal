/**
 * @file HcalMIPTracks.h
 * @brief Class used to encapsulate the results obtained from
 *        HcalMIPTracking.
 * @author Matt Solt, University of Virginia
 */

#ifndef HCAL_EVENT_HCALMIPTRACKS_H_
#define HCAL_EVENT_HCALMIPTRACKS_H_

//----------//
//   ROOT   //
//----------//
#include "TObject.h"  //For ClassDef

//----------//
//   LDMX   //
//----------//
#include "Hcal/Event/HcalHit.h"
#include "Hcal/Event/HcalMIPTrack.h"

namespace ldmx {

class HcalMIPTracks {
 public:
  /** Constructor */
  HcalMIPTracks();

  /** Destructor */
  ~HcalMIPTracks();

  /** Reset the object. */
  void Clear();

  /** Print out the object */
  void Print() const;

  inline bool getIsTriggered() const { return isTriggered_; };

  //inline std::vector<ldmx::HcalHit> getMIPTrackHits() const { return trackHits_; };

  inline int getNTracks() const { return nTracks_; };

  inline int getTriggerStart() const { return startTrig_; };

  inline int getTriggeredLayers() const { return nTrigLay_; };

  inline std::vector<ldmx::HcalMIPTrack> getMIPTracks() const { return mipTracks_; };

  /*inline float getX() const { return x_; };

  inline float getY() const { return y_; };

  inline float getDX() const { return dx_; };

  inline float getDY() const { return dy_; };

  inline float getXX() const { return xx_; };

  inline float getXY() const { return xy_; };

  inline float getXDX() const { return xdx_; };

  inline float getXDY() const { return xdy_; };

  inline float getYY() const { return yy_; };

  inline float getYDX() const { return ydx_; };

  inline float getYDY() const { return ydy_; };

  inline float getDXDX() const { return dxdx_; };

  inline float getDXDY() const { return dxdy_; };

  inline float getDYDY() const { return dydy_; };

  inline void setMIPTrackHits(const std::vector<ldmx::HcalHit> trackHits) {
    trackHits_ = trackHits;
  }*/

  inline void setNTracks(const int nTracks) {
    nTracks_ = nTracks;
  }

  inline void setTriggerStart(const int startTrig) {
    startTrig_ = startTrig;
  }

  inline void setTriggeredLayers(const int nTrigLay) {
    nTrigLay_ = nTrigLay;
  }

  /*inline void setX(const float x) {
    x_ = x;
  }

  inline void setY(const float y) {
    y_ = y;
  }

  inline void setDX(const float dx) {
    dx_ = dx;
  }

  inline void setDY(const float dy) {
    dy_ = dy;
  }

  inline void setXX(const float xx) {
    xx_ = xx;
  }

  inline void setXY(const float xy) {
    xy_ = xy;
  }

  inline void setXDX(const float xdx) {
    xdx_ = xdx;
  }

  inline void setXDY(const float xdy) {
    xdy_ = xdy;
  }

  inline void setYY(const float yy) {
    yy_ = yy;
  }

  inline void setYDX(const float ydx) {
    ydx_ = ydx;
  }

  inline void setYDY(const float ydy) {
    ydy_ = ydy;
  }

  inline void setDXDX(const float dxdx) {
    dxdx_ = dxdx;
  }

  inline void setDXDY(const float dxdy) {
    dxdy_ = dxdy;
  }

  inline void setDYDY(const float dydy) {
    dydy_ = dydy;
  }*/

  inline void setIsTriggered(const bool& isTriggered = true) {
    isTriggered_ = isTriggered;
  }

  inline void setMIPTracks(std::vector<ldmx::HcalMIPTrack> mipTracks) {
    mipTracks_ = mipTracks;
  }

 private:

  //std::vector<ldmx::HcalHit> trackHits_;

  int nTracks_{0};

  int startTrig_{0};

  int nTrigLay_{0};

  std::vector<ldmx::HcalMIPTrack> mipTracks_;

  /*float x_{-9999.};

  float y_{-9999.};

  float dx_{-9999.};

  float dy_{-9999.};

  float xx_{-9999.};

  float xy_{-9999.};

  float xdx_{-9999.};

  float xdy_{-9999.};

  float yy_{-9999.};

  float ydx_{-9999.};

  float ydy_{-9999.};

  float dxdx_{-9999.};

  float dxdy_{-9999.};

  float dydy_{-9999.};*/

  /** Flag indicating whether the event passes the Hcal veto. */
  bool isTriggered_{false};

  ClassDef(HcalMIPTracks, 2);

};  // HcalMIPTracks
}  // namespace ldmx

#endif  // HCAL_EVENT_HCALMIPTRACKS_H_
