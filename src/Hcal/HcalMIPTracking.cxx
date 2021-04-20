#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "Hcal/HcalMIPTracking.h"
#include "Hcal/HcalDigiProducer.h"
#include "Hcal/Event/HcalMIPTracks.h"
#include "Hcal/Event/HcalMIPTrack.h"

#include <exception>
#include <iostream>

/*
TODO:
Calculate Covariances correctly
Fix same hits
*/

namespace hcal {

HcalMIPTracking::HcalMIPTracking(const std::string& name,
                                   framework::Process& process)
    : Producer(name, process) {}

void HcalMIPTracking::configure(framework::config::Parameters& parameters) {
  BACK_HCAL_START_Z_ = parameters.getParameter<double>("BACK_HCAL_START_Z_");
  MIP_MIN_PE_ = parameters.getParameter<double>("MIP_MIN_PE_");
  MIP_MAX_PE_ = parameters.getParameter<double>("MIP_MAX_PE_");
  MIN_TRACK_HITS_ = parameters.getParameter<int>("MIN_TRACK_HITS_");
  MIN_SEED_HITS_ = parameters.getParameter<int>("MIN_SEED_HITS_");
  MAX_SEED_HIT_ERROR_ = parameters.getParameter<double>("MAX_SEED_HIT_ERROR_");
  MAX_TRACK_EXTRAP_SIGMA_ = parameters.getParameter<double>("MAX_TRACK_EXTRAP_SIGMA_");
  MAX_LAYERS_CONSEC_MISSED_ = parameters.getParameter<int>("MAX_LAYERS_CONSEC_MISSED_");
  USE_ISOLATED_HITS_ = parameters.getParameter<bool>("USE_ISOLATED_HITS_");
  NUM_HITS_REQ_ = parameters.getParameter<int>("NUM_HITS_REQ_");
  NUM_HITS_IN_GROUP_ = parameters.getParameter<int>("NUM_HITS_IN_GROUP_");
  NUM_GROUPS_REQ_ = parameters.getParameter<int>("NUM_GROUPS_REQ_");
  NUM_LAY_PER_GROUP_  = parameters.getParameter<int>("NUM_LAY_PER_GROUP_");
  NUM_GROUPS_PER_LAY_  = parameters.getParameter<int>("NUM_GROUPS_PER_LAY_");

  STRIPS_BACK_PER_LAYER_ = parameters.getParameter<int>("strips_back_per_layer");
  NUM_BACK_HCAL_LAYERS_ = parameters.getParameter<int>("num_back_hcal_layers");
  STRIPS_SIDE_TB_PER_LAYER_ = parameters.getParameter<int>("strips_side_tb_per_layer");
  NUM_SIDE_TB_HCAL_LAYERS_ = parameters.getParameter<int>("num_side_tb_hcal_layers");
  STRIPS_SIDE_LR_PER_LAYER_ = parameters.getParameter<int>("strips_side_lr_per_layer");
  NUM_SIDE_LR_HCAL_LAYERS_ = parameters.getParameter<int>("num_side_lr_hcal_layers");
}

void HcalMIPTracking::produce(framework::Event& event) {

    std::vector<ldmx::HcalHit> hcalHits = event.getCollection<ldmx::HcalHit>("HcalRecHits");

    std::vector<ldmx::HcalHit> hcalIsoHits = FindIsolatedHits(hcalHits, USE_ISOLATED_HITS_ );

    std::map<int, std::vector<ldmx::HcalHit>> hcalIsoSortedHits;

    //for (const ldmx::HcalHit &hit : hcalIsoHits ) {
    for (const ldmx::HcalHit &hit : hcalHits ) {
      ldmx::HcalID detID(hit.getID());
      int layer = detID.getLayerID();
      int n = layer;// - 1;
      if(hcalIsoSortedHits.count(n) < 1){
        std::vector<ldmx::HcalHit> temp;
        temp.push_back(hit);
        hcalIsoSortedHits.insert(std::pair<int, std::vector<ldmx::HcalHit>>(n, temp));
      }
      else{
        hcalIsoSortedHits[n].push_back(hit);
      }
    }

    //bool trigger = IsTriggered(hcalIsoSortedHits);

    std::vector<std::vector<ldmx::HcalHit>> tracklist = FindTracks(hcalIsoSortedHits);

    ldmx::HcalMIPTracks tracks;
    tracks.setNTracks(tracklist.size());
    std::vector<int> layers = TriggeredLayers(hcalHits);
    tracks.setTriggerStart(layers[0]);
    tracks.setTriggeredLayers(layers[1]);
    bool trigger = layers[1] >= NUM_GROUPS_REQ_;
    tracks.setIsTriggered(trigger);
    std::vector<ldmx::HcalMIPTrack> miptracks;
    for(std::vector<ldmx::HcalHit> trackhits : tracklist){
      ldmx::HcalMIPTrack track;
      track.setMIPTrackHits(trackhits);
      float *fit;
      fit = fitTrackLS(trackhits);
      track.setX(fit[0]);
      track.setY(fit[1]);
      track.setDX(fit[2]);
      track.setDY(fit[3]);
      track.setXX(fit[4]);
      track.setXY(fit[5]);
      track.setXDX(fit[6]);
      track.setXDY(fit[7]);
      track.setYY(fit[8]);
      track.setYDX(fit[9]);
      track.setYDY(fit[10]);
      track.setDXDX(fit[11]);
      track.setDXDY(fit[12]);
      track.setDYDY(fit[13]);
      miptracks.push_back(track);
    }
    tracks.setMIPTracks(miptracks);
    event.add("HcalMIPTracks", tracks);
    return;
}

//std::vector<int> HcalMIPTracking::TriggeredLayers(std::map<int, std::vector<ldmx::HcalHit>> &hitmap){
std::vector<int> HcalMIPTracking::TriggeredLayers(std::vector<ldmx::HcalHit> &hits){

  std::vector<float> temp;
  for (int i = 0; i < NUM_GROUPS_PER_LAY_; i++){
    temp.push_back(0);
  }
  std::vector<std::vector<float>> PEsum;
  for (int i = 0; i < int(NUM_BACK_HCAL_LAYERS_ / NUM_LAY_PER_GROUP_); i++){
    PEsum.push_back(temp);
  }

  for (const ldmx::HcalHit &hit : hits ) {
    ldmx::HcalID detID(hit.getID());
    int section = detID.getSection();
    if(section != 0){
      continue;
    }
    int strip = detID.getStrip();
    int layer = detID.getLayerID();
    int PE = hit.getPE();
    int block = int(layer / NUM_LAY_PER_GROUP_);
    int group = int(strip / (STRIPS_BACK_PER_LAYER_ / NUM_GROUPS_PER_LAY_));
    PEsum[block][group] = PEsum[block][group] + PE;
  }

  std::vector<bool> triggeredlayer;
  for (int i = 0; i < PEsum.size(); i++){
    bool triggered = false;
    for (int j = 0; j < PEsum[i].size(); j++){
      if(PEsum[i][j] > MIP_MIN_PE_ && PEsum[i][j] < MIP_MAX_PE_){
        triggered = true;
        break;
      }
    }
    triggeredlayer.push_back(triggered);
  }
  int nLay = 0;
  int nLaymax = 0;
  int start = -9999;
  for (int i = 0; i < triggeredlayer.size(); i++){
    if(triggeredlayer[i]){
      nLay++;
    }
    else{
      if(nLay > nLaymax){
        nLaymax = nLay;
        start = i - nLay;
      }
      nLay = 0;
    }
    if(i == triggeredlayer.size() - 1  && nLay > nLaymax){
      nLaymax = nLay;
      start = i - nLay;
    }
  }
  std::vector<int> output;
  output.push_back(start);
  output.push_back(nLaymax);
  return output;
}

std::vector<ldmx::HcalHit> HcalMIPTracking::FindIsolatedHits(std::vector<ldmx::HcalHit> &hits, bool &use_isolated){
  std::vector<ldmx::HcalHit> isohits;
  int i = 0;
  for (const ldmx::HcalHit &hit : hits ) {
    i++;
    ldmx::HcalID detID(hit.getID());
    int section = detID.getSection();
    if(section != 0){
      continue;
    }
    int strip = detID.getStrip();
    int layer = detID.getLayerID();
    int PE      = hit.getPE();
    if(PE < MIP_MIN_PE_ || PE > MIP_MAX_PE_){
      continue;
    }
    bool isolated = true;
    int j = 0;
    for (const ldmx::HcalHit &hit2 : hits ) {
      j++;
      if(i == j){ //This is the same hit
        continue;
      }
      ldmx::HcalID detID2(hit2.getID());
      int section2 = detID2.getSection();
      if(section2 != 0){
        continue;
      }
      int strip2 = detID2.getStrip();
      int layer2 = detID2.getLayerID();
      int PE2      = hit2.getPE();
      if(PE2 < MIP_MIN_PE_ || PE2 > MIP_MAX_PE_){
        continue;
      }
      if(layer == layer2 && std::abs(strip - strip2) <= 1){
        isolated = false;
        break;
      }
    }
    if(isolated || !use_isolated){
      isohits.push_back(hit);
    }
  }
  return isohits;
}

float HcalMIPTracking::CalcDist(ldmx::HcalHit &h, std::vector<ldmx::HcalHit> &hits){
  float *fit;
  fit = fitTrackLS(hits);
  double xpos = h.getXPos();
  double ypos = h.getYPos();
  double zpos = h.getZPos() - BACK_HCAL_START_Z_;
  double x0 = fit[0];
  double y0 = fit[1];
  double dx = fit[2];
  double dy = fit[3];
  double x0_err = fit[4];
  double y0_err = fit[5];
  double dx_err = fit[6];
  double dy_err = fit[7];
  double extrap_x = dx * zpos + x0;
  double extrap_y = dy * zpos + y0;
  double extrap_x_err = abs(dx) * zpos + x0_err; //I think this calc is wrong
  double extrap_y_err = abs(dy) * zpos + y0_err;
  float dist = abs(extrap_y - ypos);
  if(abs(extrap_y - ypos) / extrap_y_err > MAX_TRACK_EXTRAP_SIGMA_){
    dist = 99999;
  }
  return dist;
}

std::vector<ldmx::HcalHit> HcalMIPTracking::chooseSeed(std::vector<std::vector<ldmx::HcalHit>> &seedlist){
  std::vector<ldmx::HcalHit> bestseed = seedlist[0];
  if(seedlist.size() == 1){
    return bestseed;
  }
  else{
    float minerr = 9999;
    for(std::vector<ldmx::HcalHit> seed : seedlist){
      float *fit;
      fit = fitTrackLS(seed);
      double x0_err = fit[4];
      double y0_err = fit[5];
      double err = y0_err;
      if(err < minerr){
        bestseed = seed;
        err = minerr;
      }
    }
    return bestseed;
  }
}

std::vector<std::vector<ldmx::HcalHit>> HcalMIPTracking::FindTracks(std::map<int, std::vector<ldmx::HcalHit>> &hitmap){
  std::vector<std::vector<ldmx::HcalHit>> tracklist;
  std::vector<ldmx::HcalHit> rmhits;
  int nseeds = 0;
  for (int i = 0; i < NUM_BACK_HCAL_LAYERS_; i++){
    if(hitmap.count(i) < 1){
      continue;
    }

    std::vector<std::vector<ldmx::HcalHit>> seedlist = FindSeeds(hitmap, rmhits, i);

    nseeds = nseeds + seedlist.size();

    if(seedlist.size() < 1){
      continue;
    }
    std::vector<ldmx::HcalHit> seed = chooseSeed(seedlist);

    std::vector<ldmx::HcalHit> hitlist;

    for(const ldmx::HcalHit &hit : seed){
      hitlist.push_back(hit);
      rmhits.push_back(hit);
    }
    int seedindex = i + 4;
    int penalty = 0;
    for (int j = seedindex; j < (NUM_BACK_HCAL_LAYERS_ - 4); j++){
      if(penalty > MAX_LAYERS_CONSEC_MISSED_){
        break;
      }
      if(hitmap.count(j) < 1){
        penalty++;
        continue;
      }
      float minimum = 9999;
      float minthresh = 9998;
      ldmx::HcalHit keephit;
      for(int k = 0; k<hitmap[j].size(); k++){
         ldmx::HcalHit hit = hitmap[j][k];
         if(isOnTrack(hit, rmhits)){
           continue;
         }
         float dist = CalcDist(hit, hitlist);
         if(dist < minimum){
           minimum = dist;
           keephit = hit;
         }
      }
      if(minimum > minthresh){
        penalty++;
      }
      else{
        hitlist.push_back(keephit);
        rmhits.push_back(keephit);
        penalty = 0;
      }
    }
    tracklist.push_back(hitlist);
  }
  return tracklist;
}

std::vector<std::vector<ldmx::HcalHit>> HcalMIPTracking::FindSeeds(std::map<int, std::vector<ldmx::HcalHit>> &hitmap, std::vector<ldmx::HcalHit> &rmhit, int &index){
  std::vector<std::vector<ldmx::HcalHit>> seedlist;
  for(ldmx::HcalHit &hit1 : hitmap[index]){
    if(isOnTrack(hit1, rmhit)){
      continue;
    }
    int i = 0;
    do{
      ldmx::HcalHit hit2;
      if(hitmap[index+1].size() > 0){
        hit2 = hitmap[index+1][i];
      }
      int j = 0;
      do{
        ldmx::HcalHit hit3;
        if(hitmap[index+2].size() > 0){
          hit3 = hitmap[index+2][j];
        }
        int k = 0;
        do{
          std::vector<ldmx::HcalHit> temp;
          temp.push_back(hit1);
          if(hitmap[index+1].size() > 0 && !isOnTrack(hit2, rmhit)){
            temp.push_back(hit2);
          }
          if(hitmap[index+2].size() > 0 && !isOnTrack(hit3, rmhit)){
            temp.push_back(hit3);
          }
          ldmx::HcalHit hit4;
          if(hitmap[index+3].size() > 0 && !isOnTrack(hit3, rmhit)){
            hit4 = hitmap[index+3][k];
            temp.push_back(hit4);
          }
          seedlist.push_back(temp);
          k++;
        }
        while(hitmap[index+3].size() > k);
        j++;
      }
      while(hitmap[index+2].size() > j);
      i++;
    }
    while(hitmap[index+1].size() > i);
  }
  std::vector<std::vector<ldmx::HcalHit>> cleanseedlist = CleanSeedList(seedlist);
  return cleanseedlist;
}

bool HcalMIPTracking::isOnTrack(ldmx::HcalHit &hit, std::vector<ldmx::HcalHit> &rmhit){
  for(ldmx::HcalHit h : rmhit){
    if(abs(hit.getXPos() - h.getXPos()) < 0.001){ //This is hard-coded and needs to be fixed
      //std::cout<<std::addressof(h)<<"  "<<std::addressof(hit)<<std::endl;
      //bool issame = (h == hit);
      //std::cout<<&hit<<"  "<<&h<<std::endl;
      return true;
    }
  }
  return false;
}

std::vector<std::vector<ldmx::HcalHit>> HcalMIPTracking::CleanSeedList(std::vector<std::vector<ldmx::HcalHit>> &seedlist){
  std::vector<std::vector<ldmx::HcalHit>> cleanseedlist;
  for(int i = 0; i<seedlist.size(); i++){
    if(seedlist[i].size() < MIN_SEED_HITS_){
      continue;
    }
    float *fit;
    fit = fitTrackLS(seedlist[i]);
    double x0 = fit[0];
    double y0 = fit[1];
    double dx = fit[2];
    double dy = fit[3];
    //Update these indices
    double x0_err = fit[4];
    double y0_err = fit[5];
    double dx_err = fit[6];
    double dy_err = fit[7];
    double err = y0_err;
    std::vector<ldmx::HcalHit> badhits;
    for (const ldmx::HcalHit &hit : seedlist[i] ) {
      int PE      = hit.getPE();
      double x       = hit.getXPos();
      double y       = hit.getYPos();
      double z       = hit.getZPos();
      double extrap_x = dx * z + x0;
      double extrap_y = dy * z + y0;
      if((extrap_x - x) > MAX_SEED_HIT_ERROR_ || (extrap_y - y) > MAX_SEED_HIT_ERROR_){
        badhits.push_back(hit);
      }
    }
    if(badhits.size() > 1 || (badhits.size() > 0 && seedlist[i].size() < 4)){
      continue;
    }
    if(badhits.size() > 0){
      std::vector<ldmx::HcalHit> newlist;
      for (const ldmx::HcalHit &hit : seedlist[i] ) {
        newlist.push_back(hit);
      }
      cleanseedlist.push_back(newlist);
      continue;
    }
    cleanseedlist.push_back(seedlist[i]);
  }
  return cleanseedlist;
}

bool HcalMIPTracking::IsTriggered(std::map<int, std::vector<ldmx::HcalHit>> &hitmap){
  /*int trigger = 0;
  for (int i = 0; i < NUM_BACK_HCAL_LAYERS_; i++){
    int nHits = 0;
    if(hitmap.count(i) < 1){
      continue;
    }
    else{
      nHits++;
    }
    if(hitmap.count(i+1) > 0){
      nHits++;
    }
    if(hitmap.count(i+2) > 0){
      nHits++;
    }
    if(hitmap.count(i+3) > 0){
      nHits++;
    }
    if(nHits >= 3){
      return true;
    }
  }
  return false;*/
  std::vector<int> triggers;
  for (int i = 0; i < NUM_BACK_HCAL_LAYERS_; i++){
    int nHits = 0;
    for(int j = i; j < i+NUM_HITS_IN_GROUP_; j++){
      if(hitmap.count(j) > 0){
        nHits++;
      }
      if(nHits >= NUM_HITS_REQ_){
        triggers.push_back(1);
      }
      else{
        triggers.push_back(0);
      }
    }
  }
  for (int i = 0; i < triggers.size(); i++){
    int n = 0;
    for(int j = i; j < triggers.size(); j=j+NUM_HITS_IN_GROUP_){
      if(triggers[j] == 0){
        break;
      }
      else{
        n++;
      }
      if(n >= NUM_GROUPS_REQ_){
        return true;
      }
    }
  }
  return false;
}

float HcalMIPTracking::vectorMean(std::vector<float> &vec){
  float mean = 0;
  for(int i = 0; i < vec.size(); i++){
    mean += vec[i];
  }
  mean = mean / vec.size();
  return mean;
}

float HcalMIPTracking::vectorSD2(std::vector<float> &vec){
  float sd = 0;
  for(int i = 0; i < vec.size(); i++){
    sd += pow(vec[i], 2);
  }
  sd = (sd - vec.size() * pow(vectorMean(vec), 2)) / vec.size();
  return sd;
}

float HcalMIPTracking::vectorSS(std::vector<float> &vec1, std::vector<float> &vec2){
  float ss = 0;
  for(int i = 0; i < vec1.size(); i++){
    ss += vec1[i] * vec2[i];
  }
  ss = ss - vec1.size() * vectorMean(vec1) * vectorMean(vec2);
  return ss;
}


float* HcalMIPTracking::fitTrackLS(std::vector<ldmx::HcalHit> &hitlist){

  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;

  for (const ldmx::HcalHit &hit : hitlist ) {
    double xhit       = hit.getXPos();
    double yhit       = hit.getYPos();
    double zhit       = hit.getZPos();

    x.push_back(xhit);
    y.push_back(yhit);
    z.push_back(zhit - BACK_HCAL_START_Z_);
  }

  float x_mean = vectorMean(x);
  float y_mean = vectorMean(y);
  float z_mean = vectorMean(z);

  float sigx2 = vectorSD2(x);
  float sigy2 = vectorSD2(y);
  float sigz2 = vectorSD2(z);

  float ssxx = sigx2 * x.size();
  float ssyy = sigy2 * y.size();
  float sszz = sigz2 * z.size();

  float ssxy = vectorSS(x, y);
  float ssxz = vectorSS(x, z);
  float ssyz = vectorSS(y, z);

  float covxy = ssxy / x.size();
  float covxz = ssxz / x.size();
  float covyz = ssyz / y.size();

  float dx = covxz / sigz2;
  float x0 = x_mean - dx * z_mean;
  float sx = sqrt((ssxx - dx * ssxz) / (z.size() - 2));
  float x0_err = sx * sqrt(1 / z.size() + pow(z_mean, 2) / (sigz2 * z.size()));
  float dx_err = sx / sqrt(sigz2 * z.size());

  float dy = covyz / sigz2;
  float y0 = y_mean - dy * z_mean;
  float sy = sqrt((ssyy - dy * ssyz) / (z.size() - 2));
  float y0_err = sy * sqrt(1 / z.size() + pow(z_mean, 2) / (sigz2 * z.size()));
  float dy_err = sy / sqrt(sigz2 * z.size());

  /*float xx = vectorSS(x, x) / x.size();
  float xy = vectorSS(x, x) / x.size();
  float xdx = vectorSS(x, dx) / x.size();
  float xdy = vectorSS(x, dy) / x.size();
  float yy = vectorSS(y, y) / x.size();
  float ydx = vectorSS(y, dx) / x.size();
  float ydy = vectorSS(y, dy) / x.size();
  float dxdx = vectorSS(dx, dx) / x.size();
  float dxdy = vectorSS(dx, dy) / x.size();
  float dydy = vectorSS(dy, dy) / x.size();*/

  float * params = new float[14];
  params[0] = x0;
  params[1] = y0;
  params[2] = dx;
  params[3] = dy;
  params[4] = x0_err*x0_err;
  params[5] = 0;
  params[6] = 0.;
  params[7] = 0.;
  params[8] = y0_err*y0_err;
  params[9] = 0.;
  params[10] = 0.;
  params[11] = dx_err*dx_err;
  params[12] = 0.;
  params[13] = dy_err*dy_err;
  /*params[4] = xx;
  params[5] = xy;
  params[6] = xdx;
  params[7] = xdy;
  params[8] = yy;
  params[9] = ydx;
  params[10] = ydy;
  params[11] = dxdx;
  params[12] = dxdy;
  params[13] = dydy;*/
  return params;
}

}  // namespace hcal

DECLARE_PRODUCER_NS(hcal, HcalMIPTracking);
