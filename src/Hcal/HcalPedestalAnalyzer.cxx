/** @class Service application which creates standard CSV-format pedestal files
 */

#include "DetDescr/HcalDigiID.h"

#include <fire/Processor.h>

#include "Recon/Event/HgcrocDigiCollection.h"

namespace hcal {

class HcalPedestalAnalyzer : public fire::Processor {
  std::string input_name_, input_pass_;
  std::string output_file_, comments_;
  bool store_filtered_adcs_;
  bool filter_noTOT;
  bool filter_noTOA;
  int low_cutoff_, high_cutoff_;

  struct Channel {
    Channel() : sum{0}, sum_sq{0}, entries{0}, rejects{4, 0} {}
    /// collection of hits accumulated to produce appropriately-binned
    /// histograms
    std::vector<int> adcs;
    /// Sum of values
    uint64_t sum;
    /// Sum of values squared
    double sum_sq;
    /// Number of entries
    int entries;
    /// counts of various rejections
    std::vector<int> rejects;
  };

  std::map<ldmx::HcalDigiID, Channel> pedestal_data_;

  void create_and_fill(Channel& chan, ldmx::HcalDigiID detid);

 public:
  HcalPedestalAnalyzer(const fire::config::Parameters& ps)
    : fire::Processor(ps) {
    input_name_ = ps.get<std::string>("input_name");
    input_pass_ = ps.get<std::string>("input_pass");
    output_file_ = ps.get<std::string>("output_file");
    comments_ = ps.get<std::string>("comments");

    store_filtered_adcs_ = ps.get<bool>("store_filtered_adcs", false);
    filter_noTOT = ps.get<bool>("filter_noTOT", true);
    filter_noTOA = ps.get<bool>("filter_noTOA", true);
    low_cutoff_ = ps.get<int>("low_cutoff", 10);
    high_cutoff_ = ps.get<int>("high_cutoff", 512);
  }

  virtual void process(fire::Event& event) final override;
  virtual void onProcessEnd();
};

void HcalPedestalAnalyzer::process(fire::Event& event) {
  auto const& digis{
      event.get<ldmx::HgcrocDigiCollection>(input_name_, input_pass_)};

  for (std::size_t i_digi{0}; i_digi < digis.size(); i_digi++) {
    auto d{digis.getDigi(i_digi)};
    ldmx::HcalDigiID detid(d.id());

    Channel& chan = pedestal_data_[detid];
    chan.adcs.clear();

    bool has_tot = false;
    bool has_toa = false;
    bool has_under = false;
    bool has_over = false;

    for (int i = 0; i < digis.getNumSamplesPerDigi(); i++) {
      if (d.at(i).tot() > 0) has_tot = true;
      if (d.at(i).toa() > 0) has_toa = true;
      if (d.at(i).adc_t() < low_cutoff_) has_under = true;
      if (d.at(i).adc_t() > high_cutoff_) has_over = true;
    }

    if (has_tot && filter_noTOT) chan.rejects[0]++;
    if (has_toa && filter_noTOA) chan.rejects[1]++;
    if (has_under) chan.rejects[2]++;
    if (has_over) chan.rejects[3]++;

    if (has_tot && filter_noTOT) continue;  // ignore this
    if (has_toa && filter_noTOA) continue;  // ignore this
    if (has_under)
      continue;  // ignore this, set threshold to zero to disable requirement
    if (has_over)
      continue;  // ignore this, set threshold larger than 1024 to disable
                 // requirement

    for (int i = 0; i < digis.getNumSamplesPerDigi(); i++) {
      int adc = d.at(i).adc_t();

      chan.sum += adc;
      chan.sum_sq += adc * adc;
      chan.entries++;
      chan.adcs.push_back(adc);
    }
  }

  if (store_filtered_adcs_) {
    std::map<int,std::vector<int>> filtered_adcs;
    for (auto& [id, chan] : pedestal_data_) {
      if (chan.adcs.size() > 0) {
        filtered_adcs[id.raw()] = chan.adcs;
      }
    }
    event.add("FilteredADCs", filtered_adcs);
  }
}

void HcalPedestalAnalyzer::onProcessEnd() {
  FILE* fout = fopen(output_file_.c_str(), "w");

  time_t t = time(NULL);
  struct tm* gmtm = gmtime(&t);
  char times[1024];
  strftime(times, sizeof(times), "%Y-%m-%d %H:%M:%S GMT", gmtm);
  fprintf(fout, "# %s\n", comments_.c_str());
  fprintf(fout, "# Produced %s\n", times);
  fprintf(fout, "DetID,PEDESTAL_ADC,PEDESTAL_RMS_ADC\n");

  for (auto ichan : pedestal_data_) {
    if (ichan.second.entries == 0) {
      std::cout << "All entries filtered for " << ichan.first << " for TOT "
                << ichan.second.rejects[0] << " for TOA "
                << ichan.second.rejects[1] << " for underthreshold "
                << ichan.second.rejects[2] << " for overthreshold "
                << ichan.second.rejects[3] << std::endl;
      continue;  // all entries were filtered out
    }

    double mean = ichan.second.sum * 1.0 / ichan.second.entries;
    double rms = sqrt(ichan.second.sum_sq / ichan.second.entries - mean * mean);
    fprintf(fout, "0x%08x,%9.3f,%9.3f\n", ichan.first.raw(), mean, rms);
  }

  fclose(fout);
}

}  // namespace hcal

DECLARE_PROCESSOR(hcal::HcalPedestalAnalyzer);
