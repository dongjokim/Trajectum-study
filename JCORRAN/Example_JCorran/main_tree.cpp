#include <fstream>
#include <iostream>
#include <string>

#include <TClonesArray.h>
#include <TFile.h>
#include <TROOT.h>

#include "args.hxx"
#include "AliJBaseEventHeader.h"
#include "AliJBaseTrack.h"
#include "AliJFFlucAnalysisTProfile.h"
#include "AliAnalysisPtVn.h"
#include "AliAnalysisSPCRun2.h"
#include "JTreeDataManager.h"

int main(int argc, char **pargv) {
	args::ArgumentParser parser("JFLUC analyzer (JCORRAN Tree input)", "");
	args::ValueFlag<std::string> outfilen(parser, "dst", "Output ROOT file", {'o', "outfile"}, "AnalysisResults.root");
	args::ValueFlag<std::string> filelist(parser, "list", "Text file listing JCORRAN Tree ROOT files", {'l', "list"}, "");
	args::Positional<std::string> singleFile(parser, "input.root", "Single JCORRAN Tree file (if --list is not used)");
	args::ValueFlag<std::string> binningMode(parser, "mode", "Binning: PbPb_cent or OO_cent", {'b', "binning"}, "OO_cent");
	args::ValueFlag<double> minPt(parser, "pTmin", "Minimum track pT", {"pTmin"}, 0.2);
	args::ValueFlag<double> maxPt(parser, "pTmax", "Maximum track pT", {"pTmax"}, 3.0);
	args::ValueFlag<double> absEtaMin(parser, "absEtaMin", "Minimum |eta| for eta gap", {"absEtaMin"}, 0.4);
	args::ValueFlag<double> absEtaMax(parser, "absEtaMax", "Maximum track |eta|", {"absEtaMax"}, 0.8);
	args::ValueFlag<std::string> addObs(parser, "addObs", "Additional observables: spc", {"addObs"}, "none");

	try {
		parser.ParseCLI(argc, pargv);
	} catch (args::Help &) {
		std::cout << parser;
		return 0;
	} catch (args::ParseError &e) {
		std::cout << e.what() << "\n" << parser;
		return 1;
	}

	std::string listPath = filelist.Get();
	std::string tmpList;
	if (listPath.empty()) {
		if (!singleFile)
			return (std::cerr << "Need --list or a ROOT file.\n"), 1;
		tmpList = "/tmp/jcorran_tree_list.txt";
		std::ofstream ofs(tmpList);
		ofs << args::get(singleFile) << "\n";
		listPath = tmpList;
	}

	gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

	auto binning = AliJFFlucAnalysisTProfile::BINNING_CENT_PbPb;
	if (binningMode.Get() == "OO_cent")
		binning = AliJFFlucAnalysisTProfile::BINNING_CENT_OO;

	const bool bSPC = addObs.Get().find("spc") != std::string::npos;

	TFile *pfo = new TFile(outfilen.Get().c_str(), "recreate");
	TDirectory *phydroDir = pfo->mkdir("hydro");
	pfo->cd("hydro");

	auto *pfa = new AliJFFlucAnalysisTProfile("bayesian-hydro-jfluc");
	pfa->SetBinning(binning);
	pfa->AddFlags(AliJFFlucAnalysisTProfile::FLUC_EBE_WEIGHTING);
	pfa->SetEtaRange(absEtaMin.Get(), absEtaMax.Get());
	pfa->UserCreateOutputObjects();

	phydroDir->mkdir("vnptCorr");
	phydroDir->cd("vnptCorr");
	auto *pca_PtVn = new AliAnalysisPtVn("PtVnCorr");
	pca_PtVn->UserCreateOutputObjects();
	pca_PtVn->SetPtSubRange(absEtaMin.Get(), absEtaMax.Get());

	AliAnalysisSPCRun2 *pspc = nullptr;
	if (bSPC) {
		phydroDir->mkdir("SPC");
		phydroDir->cd("SPC");
		pspc = new AliAnalysisSPCRun2("SPC_0");
		pspc->UserCreateOutputObjects();
	}

	JTreeDataManager dm;
	dm.ChainInputStream(listPath.c_str());
	TClonesArray *pinputList = new TClonesArray("AliJBaseTrack", 2500);

	const int nEvt = dm.GetNEvents();
	int used = 0;
	for (int ie = 0; ie < nEvt; ++ie) {
		dm.LoadEvent(ie);
		auto *hdrList = dm.GetEventHeaderList();
		if (!hdrList || hdrList->GetEntriesFast() < 1)
			continue;
		auto *hdr = static_cast<AliJBaseEventHeader *>(hdrList->At(0));
		const float cent = hdr->GetCentrality();
		const int fCBin = AliJFFlucAnalysisTProfile::GetBin(cent, binning);
		if (fCBin < 0)
			continue;

		pinputList->Clear();
		dm.RegisterList(pinputList, nullptr);

		TClonesArray filtered("AliJBaseTrack", 2500);
		int nkeep = 0;
		for (int i = 0; i < pinputList->GetEntriesFast(); ++i) {
			auto *tr = static_cast<AliJBaseTrack *>(pinputList->At(i));
			const double pt = tr->Pt();
			const double eta = tr->Eta();
			if (pt < minPt.Get() || pt > maxPt.Get() || std::fabs(eta) > absEtaMax.Get())
				continue;
			new (filtered[nkeep++]) AliJBaseTrack(*tr);
		}

		pfa->SetInputList(&filtered);
		pfa->SetEventCentrality(cent);
		pfa->SetEventImpactParameter(-1);
		pfa->UserExec("");

		phydroDir->cd("vnptCorr");
		pca_PtVn->SetInputList(&filtered);
		pca_PtVn->SetEventCentrality(cent);
		pca_PtVn->UserExec("");

		if (bSPC) {
			phydroDir->cd("SPC");
			pspc->SetInputList(&filtered);
			pspc->SetEventCentrality(cent);
			pspc->SetSPC(0);
			pspc->UserExec("");
		}
		++used;
	}

	phydroDir->cd("vnptCorr");
	pca_PtVn->WriteLists("OutputPtVnCorr");
	if (bSPC) {
		phydroDir->cd("SPC");
		pspc->WriteLists("OutputListSPC");
	}
	pfo->Write();
	pfo->Close();
	std::cout << "analysis DONE: " << used << " / " << nEvt << " events -> " << outfilen.Get() << "\n";
	return 0;
}
