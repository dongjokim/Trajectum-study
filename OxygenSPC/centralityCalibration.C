#include "src/AliJBaseEventHeader.h"
#include "src/AliJBaseTrack.h"
#include <cmath>

R__LOAD_LIBRARY(src/AliJBaseEventHeader_cxx.so);
R__LOAD_LIBRARY(src/AliJBaseTrack_cxx.so);

void centralityCalibration() {
  TString infile_path = "../OxygenHydro-Sim/rootfiles/WS/merged";
  TChain *jTree = new TChain("jTree");
  jTree->Add(infile_path + "/*.root");

  TClonesArray *event = new TClonesArray("AliJBaseEventHeader");
  TClonesArray *tracks = new TClonesArray("AliJBaseTrack");

  jTree->SetBranchAddress("JEventHeaderList", &event);
  jTree->SetBranchAddress("JTrackList", &tracks);

  int nevents = jTree->GetEntries();
  cout << "Total events found = " << nevents << endl;

  TH1D *h = new TH1D("nChFT0M","",3000,0,3000);

  for(int i=0; i<nevents; i++){
	if(!(nevents<1000)){
	if(i%10 == 0) { int value = TMath::Log10(i)+1; for(int cc=1; cc<=value; cc++){cout<<"\b";} cout<<i<<flush;}
	} 
    event->Clear();
    tracks->Clear();
    jTree->GetEntry(i);

    AliJBaseEventHeader *hdr = (AliJBaseEventHeader*)event->At(0);
    int nChFT0M = hdr->GetYVertexMC();
    // int nChFT0M = hdr->GetZertexErr();
    h->Fill(nChFT0M);
  }

  h->Draw();

  TFile *outfile = new TFile("centralityCalibrationHist.root", "recreate");
  outfile->cd();
  h->Write();
  outfile->Close();

  vector<double> cent_bins = {0.0, 30, 40, 50, 60, 70, 80, 90, 95, 100};
  int cent_bin_size = std::size(cent_bins)-1;
  cout<<"Running calibration ..."<<endl;

  double total = h->Integral();
  cout<<total<<endl;

  for(int i=1; i<= cent_bin_size; i++){
    double current_integral = 0;
    double ratio = 0;
    for(int bin=1; bin<=h->GetNbinsX(); bin++){
        current_integral += h->GetBinContent(bin);
        ratio = (current_integral / total) * 100.0;
        if(ratio > cent_bins[i]){
            double bin_center = h->GetBinCenter(bin-1);
            cout << "Centrality bin " << cent_bins[i-1]<<" "<< cent_bins[i] << ": Nch cut = " << bin_center 
                 << ", Integral = " << current_integral 
                 << ", Ratio = " << ratio << "%" << endl;
            break;
        }
    }
  }  
  

}
