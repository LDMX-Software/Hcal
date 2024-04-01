#include "Hcal/HcalResolutionAnalyzer.h"

namespace hcal {

void HcalResolutionAnalyzer::analyze(const framework::Event& event) {
  const auto& hcalrecHits{
      event.getCollection<ldmx::HcalHit>(rec_coll_name, input_pass_name)};
  const auto& hcalSimHits{event.getCollection<ldmx::SimCalorimeterHit>(
      sim_coll_name, input_pass_name)};
  double total_energy{0.};
  for (const auto& hit : hcalSimHits) {
    auto e{hit.getEdep()};
    total_energy += e;
  }
  histogram_sim->Fill(total_energy);
  total_energy = 0;
  for (const auto& hit : hcalrecHits) {
    auto e{hit.getEnergy()};
    total_energy += e;
  }
  histogram_rec->Fill(total_energy);
  if (events_processed_this_run++ > max_events_per_run) {
    abortEvent();
  }
}

void HcalResolutionAnalyzer::onFileClose(framework::EventFile& eventFile) {
  run_results_sim[current_particle].push_back(
      HcalResolutionAnalyzer::RunResult{histogram_sim, current_energy});
  run_results_rec[current_particle].push_back(
      HcalResolutionAnalyzer::RunResult{histogram_rec, current_energy});
}

void HcalResolutionAnalyzer::onFileOpen(framework::EventFile& eventFile) {
  events_processed_this_run = 0;
  auto runHeader = eventFile.getRunHeader(run_number);
  current_energy = getCurrentBeamEnergy(runHeader);
  current_particle = getCurrentParticle(runHeader);
  auto [sim_key, rec_key] = getHistogramKeys(current_particle, current_energy);
  histograms_.create(sim_key, "Total energy deposited [MeV]", bins, 0,
                     current_energy * max_energy_histogram_fraction);
  histograms_.create(rec_key, "Total reconstructed energy [MeV]", bins, 0,
                     current_energy * max_energy_histogram_fraction);

  histogram_sim = histograms_.get(sim_key);
  histogram_rec = histograms_.get(rec_key);
}

void HcalResolutionAnalyzer::onProcessEnd() {
  getHistoDirectory();

  // For each particle type that we have looked at
  // - Create a graph for resolution as function of beam energy
  // - Create a graph for linearity as function of beam energy
  // - Create a multigraph for resolution with all particles at once
  // - Create a multigraph for linearity with all particles at once

  // Find all particles that we have looked at from the maps
  auto particles{getProcessedParticles()};
  auto energies{getProcessedEnergies()};
  const auto N{static_cast<int>(energies.size())};

  int marker_counter{20};
  auto fillGraphs =
      [&](const std::vector<HcalResolutionAnalyzer::RunResult>& results) {
        std::vector<double> beam_energy(N);
        std::vector<double> beam_energy_err(N);
        std::vector<double> linearity(N);
        std::vector<double> linearity_err(N);
        std::vector<double> resolution(N);
        std::vector<double> resolution_err(N);
        for (int i{0}; i < static_cast<int>(results.size()); ++i) {
          auto result = results[i];
          beam_energy[i] = result.beam_energy;
          beam_energy_err[i] = 0;
          linearity[i] = result.mean;
          linearity_err[i] = result.mean_err;
          resolution[i] = result.resolution;
          resolution_err[i] = result.resolution_err;
        }
        auto lin =
            new TGraphErrors(N, beam_energy.data(), linearity.data(),
                             beam_energy_err.data(), linearity_err.data());
        auto res =
            new TGraphErrors(N, beam_energy.data(), resolution.data(),
                             beam_energy_err.data(), resolution_err.data());
        lin->GetXaxis()->SetTitle("E_{beam} [MeV]");
        lin->GetYaxis()->SetTitle("\\mu_E [MeV]");
        lin->SetLineStyle(0);
        lin->SetMarkerStyle(marker_counter++);
        res->GetXaxis()->SetTitle("E_{beam} [MeV]");
        res->GetYaxis()->SetTitle("\\frac{\\sigma(E)}{\\mu_E}");
        res->SetLineStyle(0);
        res->SetMarkerStyle(marker_counter++);
        // Wish that TGraph had a nicer interface....
        return std::make_pair(lin, res);
      };
  TMultiGraph* sim_multigraph_resolution = new TMultiGraph();
  TMultiGraph* sim_multigraph_linearity = new TMultiGraph();
  TMultiGraph* rec_multigraph_resolution = new TMultiGraph();
  TMultiGraph* rec_multigraph_linearity = new TMultiGraph();

  TGraphErrors* sim_e{nullptr};
  TGraphErrors* sim_h{nullptr};
  TGraphErrors* rec_e{nullptr};
  TGraphErrors* rec_h{nullptr};
  for (auto particle : particles) {
    auto sim_results{run_results_sim[particle]};
    auto rec_results{run_results_rec[particle]};
    auto [sim_linearity, sim_resolution] = fillGraphs(sim_results);
    auto [rec_linearity, rec_resolution] = fillGraphs(rec_results);
    sim_linearity->SetNameTitle(
        ("hcal_sim_linearity_" + particle).c_str(),
        ("Hcal energy response for " + particle).c_str());
    sim_resolution->SetNameTitle(
        ("hcal_sim_resolution_" + particle).c_str(),
        ("Hcal energy resolution for " + particle).c_str());
    rec_linearity->SetNameTitle(
        ("hcal_rec_linearity_" + particle).c_str(),
        ("Hcal energy response for " + particle).c_str());
    rec_resolution->SetNameTitle(
        ("hcal_rec_resolution_" + particle).c_str(),
        ("Hcal energy resolution for " + particle).c_str());
    sim_multigraph_resolution->Add(sim_resolution);
    sim_multigraph_linearity->Add(sim_linearity);
    rec_multigraph_resolution->Add(rec_resolution);
    rec_multigraph_linearity->Add(rec_linearity);

    if (particle == "e-") {
      sim_e = sim_linearity;
      rec_e = rec_linearity;
    } else if (particle == "pi-") {
      sim_h = sim_linearity;
      rec_h = rec_linearity;
    } else if (particle == "pi+") {
      sim_h = sim_linearity;
      rec_h = rec_linearity;
    }

    sim_linearity->Write();
    sim_resolution->Write();
    rec_linearity->Write();
    rec_resolution->Write();
  }
  sim_multigraph_resolution->SetNameTitle(
      "hcal_sim_resolution", "Hcal energy resolution for all particles");
  sim_multigraph_linearity->SetNameTitle(
      "hcal_sim_linearity", "Hcal energy response for all particles");
  rec_multigraph_resolution->SetNameTitle(
      "hcal_rec_resolution", "Hcal energy resolution for all particles");
  rec_multigraph_linearity->SetNameTitle(
      "hcal_rec_linearity", "Hcal energy response for all particles");
  sim_multigraph_resolution->Write();
  sim_multigraph_linearity->Write();
  rec_multigraph_resolution->Write();
  rec_multigraph_linearity->Write();
  if (sim_e != nullptr && sim_h != nullptr) {
    auto sim_eoverh = calculateEOverH(sim_e, sim_h);
    auto rec_eoverh = calculateEOverH(rec_e, rec_h);
    sim_eoverh->SetNameTitle("hcal_sim_eoverh", "Hcal e/h");
    rec_eoverh->SetNameTitle("hcal_rec_eoverh", "Hcal e/h");
    sim_eoverh->Write();
    rec_eoverh->Write();
  }
}

TGraphErrors* HcalResolutionAnalyzer::calculateEOverH(TGraphErrors* e,
                                                      TGraphErrors* h) const {
  int N = e->GetN();

  std::vector<double> eX(N), eY(N), eXerrors(N), eYerrors(N);

  for (int i = 0; i < N; ++i) {
    // I hate this so much
    double xE, yE, xH, yH;
    e->GetPoint(i, xE, yE);
    h->GetPoint(i, xH, yH);
    // Why does TGraphErrors have a useful interface but TGraph doesn't...
    double yEerror = e->GetErrorY(i);
    double yHerror = h->GetErrorY(i);

    // Assuming that beam energies are the same because doing otherwise would be
    // wild
    eX[i] = xE;
    eXerrors[i] = 0;
    if (yH != 0) {
      eY[i] = yE / yH;
      // Error propagation, assuming (hilariously incorrectly) that e and h are
      // not correlated
      eYerrors[i] = eY[i] * std::sqrt(std::pow(yEerror / yE, 2) +
                                      std::pow(yHerror / yH, 2));
    } else {
      eY[i] = 0;
      eYerrors[i] = 0;
    }
  }

  auto eoverh{new TGraphErrors{N, eX.data(), eY.data(), eXerrors.data(),
                               eYerrors.data()}};
  eoverh->SetNameTitle("eoverh", "e/h");
  eoverh->GetXaxis()->SetTitle("E_{beam} [MeV]");
  eoverh->GetYaxis()->SetTitle("e/h");
  eoverh->SetMarkerStyle(25);
  eoverh->SetMarkerColor(kBlack);
  eoverh->SetLineStyle(0);
  return eoverh;
}
void HcalResolutionAnalyzer::configure(framework::config::Parameters& ps) {
  run_number = ps.getParameter<int>("run_number");
  sim_coll_name = ps.getParameter<std::string>("sim_coll_name");
  rec_coll_name = ps.getParameter<std::string>("rec_coll_name");
  input_pass_name = ps.getParameter<std::string>("input_pass_name");
  max_events_per_run = ps.getParameter<int>("max_events_per_run");
  max_energy_histogram_fraction =
      ps.getParameter<double>("max_energy_histogram_fraction");
  bins = ps.getParameter<int>("bins");
}

}  // namespace hcal
DECLARE_ANALYZER_NS(hcal, HcalResolutionAnalyzer);
