#ifndef HCALRESOLUTIONANALYZER_H_
#define HCALRESOLUTIONANALYZER_H_
#include <iomanip>
#include <iostream>
#include <sstream>

#include "DetDescr/HcalGeometry.h"
#include "DetDescr/HcalID.h"
#include "Framework/EventFile.h"
#include "Framework/EventProcessor.h"
#include "Hcal/Event/HcalHit.h"
#include "SimCore/Event/SimCalorimeterHit.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
namespace hcal {
class HcalResolutionAnalyzer : public framework::Analyzer {
 public:
  HcalResolutionAnalyzer(const std::string& name, framework::Process& process)
      : framework::Analyzer{name, process} {}
  void onProcessEnd() override;
  virtual ~HcalResolutionAnalyzer() = default;

  std::pair<std::string, std::string> getHistogramKeys(
      const std::string& particle, const double energy) const {
    auto key{particle + "_" + std::to_string(static_cast<int>(energy))};
    return {key, "reco_" + key};
  }
  std::map<std::string, TGraphErrors*> linearityGraphs{};
  std::map<std::string, TGraph*> resolutionGraphs{};
  TH1* histogram_sim{nullptr};
  TH1* histogram_rec{nullptr};
  struct RunResult {
    double beam_energy{};
    double mean{};
    double mean_err{};
    double std{};
    double std_err{};
    double resolution{};
    double resolution_err{};
    RunResult(TH1* histogram, double energy)
        : beam_energy{energy},
          mean{histogram->GetMean()},
          mean_err{histogram->GetMeanError()},
          // Note: RMS and StdDev mean the same thing in ROOT
          std{histogram->GetStdDev()},
          std_err{histogram->GetStdDevError()} {
      resolution = std / mean;
      resolution_err = resolution * std::sqrt(std::pow(std_err / std, 2) +
                                              std::pow(mean_err / mean, 2));
    }
    bool operator<(const RunResult& other) const {
      return beam_energy < other.beam_energy;
    }
  };
  std::vector<std::string> getProcessedParticles() const {
    std::vector<std::string> particles{};
    for (const auto& [key, _] : run_results_sim) {
      particles.push_back(key);
    }
    return particles;
  }
  std::vector<double> getProcessedEnergies() const {
    // Find all energies that we have looked at from the maps
    auto temporary_result = (*run_results_sim.begin()).second;
    // We assume that each particle species has the same set of energies
    std::vector<double> energies = {};
    for (auto result : temporary_result) {
      energies.push_back(result.beam_energy);
    }
    return energies;
  }

  std::map<std::string, std::vector<HcalResolutionAnalyzer::RunResult>>
      run_results_sim{};
  std::map<std::string, std::vector<HcalResolutionAnalyzer::RunResult>>
      run_results_rec{};
  void configure(framework::config::Parameters& ps) override;
  void analyze(const framework::Event& event) override;
  void onFileOpen(framework::EventFile& eventFile) override;
  void onFileClose(framework::EventFile& eventFile) override;
  std::string getCurrentParticle(const ldmx::RunHeader& rh) const {
    return rh.getStringParameter("Gen0 Particle");
  }
  double getCurrentBeamEnergy(const ldmx::RunHeader& rh) const {
    // Convert to MeV
    return rh.getFloatParameter("Gen0 Energy [GeV]") * 1000;
  }

  TGraphErrors* calculateEOverH(TGraphErrors* e, TGraphErrors* h) const;

 private:
  int events_processed_this_run{0};
  int max_events_per_run{1000};
  int bins{100};
  double max_energy_histogram_fraction{};
  int run_number{1};
  std::string rec_coll_name{"HcalRecHits"};
  std::string sim_coll_name{"HcalRecHits"};
  std::string input_pass_name{""};
  double current_energy;
  std::string current_particle;
};
}  // namespace hcal

#endif  // HCALRESOLUTIONANALYZER_H_
