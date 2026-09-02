#include <cmath>
//root includes
#include <TFile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TProfile.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
// O2 Framework includes
#include <Framework/HistogramSpec.h>
#include <Framework/StringHelpers.h>
// Common O2 includes
#include <PWGCF/DataModel/CorrelationsDerived.h>
#include <PWGCF/JCorran/DataModel/JCatalyst.h>
#include <PWGCF/JCorran/Core/FlowJSPCAnalysis.h>
#include <PWGCF/JCorran/Core/FlowJSPCObservables.h>
// Local includes
#include "../src/AliJBaseEventHeader.h"
#include "../src/AliJBaseTrack.h"

using namespace o2;
using namespace o2::framework;
using std::cout;
using std::cerr;
using std::endl;

// Adapter so TClonesArray of AliJBaseTrack works with O2 JQVectors::Calculate
// (expects range-iterable container and track.eta() / track.phi())
struct AliJTrackAdapter {
  AliJBaseTrack* trk = nullptr;
  float eta() const { return trk ? static_cast<float>(trk->Eta()) : 0.f; }
  float phi() const { return trk ? static_cast<float>(trk->Phi()) : 0.f; }
};

TH1F* hpT = new TH1F("pt", "", 100, 0, 10);
TH1F* hphi = new TH1F("phi", "", 100, 0.0, 2 * TMath::Pi());
TH1F* heta = new TH1F("eta", "", 100, -3.0, 3.0);
TH1F* hmult = new TH1F("multiplicity", "", 100, 0.0, 500);
TH1F* hEventCounts = new TH1F("hEventCounts", "Event Statistics;Event Type;Count", 3, 0, 3);

// Single pass: filter tracks, fill QA (pT/phi/eta + SPC phi), build adapter pointing at trackList (no clone).
static void processTracksAndFillSPCQA(
  TClonesArray const* trackList,
  std::vector<AliJTrackAdapter>& trackAdapter,
  const float cfgPtMin, const float cfgPtMax, const float cfgEtaMax,
  const int cBin, const double wNUA,
  FlowJSPCAnalysis& spcAnalysis)
{
  trackAdapter.clear();
  if (!trackList) return;
  const Int_t n = trackList->GetEntriesFast();
  trackAdapter.reserve(static_cast<size_t>(n));
  for (Int_t i = 0; i < n; i++) {
    AliJBaseTrack* track = static_cast<AliJBaseTrack*>(trackList->At(i));
    if (!track) continue;
    if (std::abs(track->Eta()) > cfgEtaMax || track->Pt() < cfgPtMin || track->Pt() > cfgPtMax) continue;
    hpT->Fill(track->Pt());
    double phi = track->Phi();
    if (phi < 0) phi += 2. * TMath::Pi();
    hphi->Fill(phi);
    heta->Fill(track->Eta());
    spcAnalysis.fillQAHistograms(cBin, phi, wNUA);
    trackAdapter.push_back({track});
  }
}

int FindMultiplicityBin(int multiplicity, const std::vector<double> &binEdges) {
  for (size_t i = 0; i < binEdges.size() - 1; ++i) {
    if (multiplicity >= binEdges[i] && multiplicity < binEdges[i + 1]) {
      return i; 
    }
  }
  return -1; 
}

// Same mid-rapidity centrality as JCORRAN: |eta|<=0.5 mapped via OO5360_main_cents.csv.
// SPC histogram slots: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-100.
int FindCentralityBin(float cent) {
  const double edges[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 100};
  for (int i = 0; i < 9; ++i) {
    if (cent >= edges[i] && cent < edges[i + 1])
      return i;
  }
  if (cent >= 100.f)
    return 8;
  return -1;
}

// Runs Q-vector calculation and correlators only (tracks already processed in single loop).
// If tQVecSec / tCorrSec are non-null, elapsed seconds are accumulated for each stage.
static void runSPCCorrelators(std::vector<AliJTrackAdapter> const& trackAdapter,
  const int cfgMultMin, const int cBin, const float cfgEtaMax,
  FlowJSPCAnalysis& spcAnalysis, FlowJSPCAnalysis::JQVectorsT& jqvecs,
  double* tQVecSec = nullptr, double* tCorrSec = nullptr)
{
  if (cBin < 0 || cBin > 8) return;
  if (static_cast<int>(trackAdapter.size()) < cfgMultMin) return;
  TStopwatch sw;
  if (tQVecSec) sw.Start();
  jqvecs.Calculate(trackAdapter, 0.0f, cfgEtaMax);
  if (tQVecSec) {
    sw.Stop();
    *tQVecSec += sw.RealTime();
  }
  if (tCorrSec) sw.Start();
  spcAnalysis.setQvectors(&jqvecs);
  spcAnalysis.calculateCorrelators(cBin);
  if (tCorrSec) {
    sw.Stop();
    *tCorrSec += sw.RealTime();
  }
}

void JCorrSPCRun3(TString inputfile = "input_trees.txt",
                       TString outputfile = "correlation_hist",
                      const int cfgWhichSPC = 0, const int cfgMultMin = 6
                  ) {

  const float cfgPtMin = 0.2f, cfgPtMax = 5.0f, cfgEtaMax = 0.8f;
  const float cfgZvtxMax = 10.0f;

  // Ascending multiplicity edges for FindMultiplicityBin (same thresholds from OxygenPythia)
  const std::vector<double> multEdges = {0, 25, 40, 58, 83, 113, 152, 202, 238, 3000}; //OxygenPythia
  // const std::vector<double> multEdges = {0, 76.5, 117.5, 166.5, 228.5, 307.5, 403.5, 547.5, 661.5, 3000}; //AC-Hydro
  const std::map<int, float> cBinMap = {{0,2.5},{1,7.5},{2,15.0},{3,25.0},{4,35.0},{5,45.0},{6,55.0},{7,65.0},{8,85.0}};
  const double wNUA = 1.0;
  std::vector<AliJTrackAdapter> trackAdapter;
  trackAdapter.reserve(4096);

  cout<<"Running with cfgWhichSPC = "<<cfgWhichSPC<<
  ",cfgPtMin = "<<cfgPtMin<<", cfgPtMax = "<<cfgPtMax<<", cfgEtaMax = "<<cfgEtaMax<<", cfgZvtxMax = "<<cfgZvtxMax<<", cfgMultMin = "<<cfgMultMin<<endl;
  cout<<"Input file list at: "<<inputfile<<", Output file: "<<outputfile<<endl;

  TStopwatch totalWall;
  totalWall.Start();

  double tSpcInit = 0, tChainSetup = 0;
  double tGetEntry = 0, tProcessTracks = 0, tQVec = 0, tSPCCorr = 0;
  double tWriteOut = 0;
  TStopwatch sw;

  sw.Start();
  HistogramRegistry spcHistograms{"SPCResults", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  FlowJSPCAnalysis spcAnalysis;
  FlowJSPCAnalysis::JQVectorsT jqvecs;
  FlowJSPCObservables spcObservables;

  //init spc analysis
  spcAnalysis.setHistRegistry(&spcHistograms);
  spcAnalysis.createHistos();

  spcObservables.setSPCObservables(cfgWhichSPC);
  spcAnalysis.setFullCorrSet(spcObservables.harmonicArray);
  sw.Stop();
  tSpcInit = sw.RealTime();

  sw.Start();
  // TString infile_path = "../rootfiles";
  TChain *jTree = new TChain("jTree");
  Long64_t nFiles = jTree->Add(inputfile+"/*.root");
  if(nFiles == 0) {
    cerr << "Error: No files found matching pattern " << inputfile << "/*.root" << endl;
    totalWall.Stop();
    delete jTree;
    return;
  }

  TClonesArray *event = new TClonesArray("AliJBaseEventHeader");
  TClonesArray *trackList = new TClonesArray("AliJBaseTrack");

  jTree->SetBranchAddress("JEventHeaderList", &event);
  jTree->SetBranchAddress("JTrackList", &trackList);

  Long64_t numberEvents = jTree->GetEntries() > 15000000 ? 15000000 : jTree->GetEntries();
  cout<<"Total number of events: "<<numberEvents<<endl;
  int ieout = numberEvents / 20 > 1 ? numberEvents / 20 : 1;
  int EventCounter = 0;
  int RejectedEventCounter = 0;
  int centCounts[9] = {0};  // events per centrality bin (for debug)

  // Set bin labels for event counts histogram
  hEventCounts->GetXaxis()->SetBinLabel(1, "All Events");
  hEventCounts->GetXaxis()->SetBinLabel(2, "Rejected Events");
  hEventCounts->GetXaxis()->SetBinLabel(3, "Analyzed Events");
  sw.Stop();
  tChainSetup = sw.RealTime();

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Loop over all events
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TStopwatch wallLoop;
  wallLoop.Start();

  for (Long64_t evt = 0; evt < numberEvents; evt++) {
    if (evt % ieout == 0) {
      cout << "Processing: " << int(float(evt) / numberEvents * 100) << "%"
           << endl;
    }

    // Count all events
    hEventCounts->Fill(0.5); // Bin 1: All Events (center of bin 1)

    event->Clear();
    trackList->Clear();
    sw.Start();
    jTree->GetEntry(evt);
    sw.Stop();
    tGetEntry += sw.RealTime();

    AliJBaseEventHeader *hdr = (AliJBaseEventHeader *)event->At(0);
    if(!hdr) {
      cerr<<"Event header not found at event "<<evt<<endl;
      hEventCounts->Fill(1.5); // Bin 2: Rejected Events
      RejectedEventCounter++;
      continue;
    }
    const int cBin = FindCentralityBin(hdr->GetCentrality());

    if (cBin < 0 || cBin > 8) {
        cerr << "Centrality bin out of range" << endl;
        hEventCounts->Fill(1.5); // Bin 2: Rejected Events
        RejectedEventCounter++;
        continue;
    }

    spcHistograms.fill(HIST("FullCentrality"), cBinMap.at(cBin));
    if (cBin >= 0 && cBin < 9) centCounts[cBin]++;
    sw.Start();
    processTracksAndFillSPCQA(trackList, trackAdapter, cfgPtMin, cfgPtMax, cfgEtaMax, cBin, wNUA, spcAnalysis);
    sw.Stop();
    tProcessTracks += sw.RealTime();

    hmult->Fill(hdr->GetXVertexMC());
    runSPCCorrelators(trackAdapter, cfgMultMin, cBin, cfgEtaMax, spcAnalysis, jqvecs, &tQVec, &tSPCCorr);

    EventCounter++;
    hEventCounts->Fill(2.5); // Bin 3: Analyzed Events (center of bin 3)
    
    trackList->Clear();
  } //event loop ends
  wallLoop.Stop();

  // Debug: event count per centrality (Centrality_8 often has few events -> empty/fragile TProfile)
  cout << "[DEBUG] Events per centrality bin:";
  for (int c = 0; c < 9; c++) cout << " C" << c << "=" << centCounts[c];
  cout << endl;

  // save and exit
  sw.Start();
  outputfile += ".root";
  TFile *pfo = new TFile(outputfile.Data(), "recreate");
  // Structure: SPCResults/Centrality_0/fResults, fCovResults, phiBefore, phiAfter (and 1..8)
  pfo->mkdir("SPCResults");
  pfo->cd("SPCResults");
  if (auto h = spcHistograms.get<TH1>(HIST("FullCentrality")); h) h->Write();
  for (int c = 0; c < 9; c++) {
    TString centDir = TString::Format("Centrality_%d", c);
    if (!gDirectory->Get(centDir.Data())) gDirectory->mkdir(centDir.Data());
    gDirectory->cd(centDir.Data());

    // Debug Centrality_8: check fResults/fCovResults for empty or NaN/Inf (causes Inf/NaN in TCanvas)
    // if (c == 8) {
    //   auto pr = spcHistograms.get<TProfile>(HIST("Centrality_8/fResults"));
    //   auto pc = spcHistograms.get<TProfile>(HIST("Centrality_8/fCovResults"));
    //   auto checkProfile = [](TProfile* p, const char* name) {
    //     if (!p) return;
    //     int nb = p->GetNbinsX();
    //     int nEmpty = 0, nNan = 0, nInf = 0;
    //     for (int b = 1; b <= nb; b++) {
    //       double v = p->GetBinContent(b);
    //       if (p->GetBinEntries(b) == 0) nEmpty++;
    //       if (std::isnan(v)) nNan++;
    //       if (std::isinf(v)) nInf++;
    //     }
    //     cout << "[DEBUG] Centrality_8/" << name << ": entries=" << p->GetEntries()
    //          << " bins=" << nb << " emptyBins=" << nEmpty << " NaN=" << nNan << " Inf=" << nInf << endl;
    //   };
    //   checkProfile(pr.get(), "fResults");
    //   checkProfile(pc.get(), "fCovResults");
    // }

    // Write profiles; then phiBefore/phiAfter.
    switch (c) {
      case 0: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_0/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_0/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_0/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_0/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 1: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_1/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_1/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_1/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_1/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 2: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_2/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_2/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_2/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_2/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 3: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_3/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_3/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_3/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_3/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 4: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_4/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_4/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_4/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_4/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 5: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_5/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_5/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_5/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_5/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 6: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_6/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_6/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_6/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_6/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 7: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_7/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_7/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_7/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_7/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
      case 8: {
        auto pr = spcHistograms.get<TProfile>(HIST("Centrality_8/fResults"));
        auto pc = spcHistograms.get<TProfile>(HIST("Centrality_8/fCovResults"));
        if (pr) pr->Write("fResults");
        if (pc) pc->Write("fCovResults");
        if (auto pb = spcHistograms.get<TH1>(HIST("Centrality_8/phiBefore")); pb) pb->Write("phiBefore");
        if (auto pa = spcHistograms.get<TH1>(HIST("Centrality_8/phiAfter")); pa) pa->Write("phiAfter");
        break;
      }
    }
    pfo->cd("SPCResults");
  }
  pfo->cd();  // back to file root
  pfo->mkdir("QA");
  pfo->cd("QA");
  hEventCounts->Write();
  hpT->Write();
  hphi->Write();
  heta->Write();
  hmult->Write();

  // outputfile->Close();

  sw.Stop();
  tWriteOut = sw.RealTime();
  totalWall.Stop();
  const double tTotal = totalWall.RealTime();
  const double tLoopWall = wallLoop.RealTime();
  const double tLoopTracked = tGetEntry + tProcessTracks + tQVec + tSPCCorr;
  const double tLoopOverhead = tLoopWall - tLoopTracked;

  cout << "\n=== Timing breakdown (real time, seconds) ===\n";
  cout << std::fixed << std::setprecision(4);
  cout << "  SPC init (registry, histos, correlator set): " << tSpcInit << "\n";
  cout << "  Chain + branches + event count:             " << tChainSetup << "\n";
  cout << "  Event loop (wall):                          " << tLoopWall << "\n";
  cout << "    - TTree::GetEntry (sum):                    " << tGetEntry << "\n";
  cout << "    - processTracksAndFillSPCQA (sum):          " << tProcessTracks << "\n";
  cout << "    - JQVectors::Calculate (sum):               " << tQVec << "\n";
  cout << "    - setQvectors + calculateCorrelators (sum): " << tSPCCorr << "\n";
  cout << "    - Sum of loop stages above:                 " << tLoopTracked << "\n";
  cout << "    - Rest (selection, hist fills, overhead):   " << tLoopOverhead << "\n";
  cout << "  Output write (ROOT file + QA histos):         " << tWriteOut << "\n";
  cout << "  Total wall (JCorrSPCRun3):                    " << tTotal << "\n";
  cout << std::defaultfloat;

  cout << "Total Number of Event scanned from input = " << numberEvents << endl;
  cout << "Total Number of Event rejected           = " << RejectedEventCounter << endl;
  cout << "Total Number of Event used for analysis  = " << EventCounter << endl;
  cout << "All files properly closed. Good Bye!" << endl;
  pfo->Close();
  delete pfo;
  delete trackList;
  delete event;
  delete jTree;
}


int main(int argc, char** argv) {
    TString inputfile = "input_trees.txt";
    TString outputfile = "corr-hist";

    int cfgMultMin = 6; // min mult 6 for 3SPC, 8 for 4SPC
    int cfgWhichSPC = 0; 

    // Allow command line arguments to override defaults
    if (argc > 1) inputfile = argv[1];
    if (argc > 2) outputfile = argv[2];
    if (argc > 3) cfgWhichSPC = atoi(argv[3]);
    if (argc > 4) cfgMultMin = atoi(argv[4]);
    
    JCorrSPCRun3(inputfile, outputfile, cfgWhichSPC, cfgMultMin);

    return 0;
}
