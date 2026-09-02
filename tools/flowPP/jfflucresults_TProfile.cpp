
#include "observables.h"
#include <functional>

typedef unsigned int uint;

//#define NH 10
#define NH 13
#define NHH 6
#define NK 5
#define NCWDE 2 //wide bins

#define OO true

#if OO
	#define NCNAT 12 //native OO
	#define NC 11 //rebinned OO
	static const double CentBinsNat[NCNAT+1] = {0, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50}; //native OO
	static const double CentBins[NC+1] = {0, 2, 5 ,10, 15, 20, 25, 30, 35, 40, 45, 50}; //rebinned OO
	static const uint CentRebin[NC+1] = {0,2,3,4,5,6,7,8,9,10,11,12}; //native -> rebinned mapping
#else
	#define NCNAT 9 //native PbPb
	#define NC 7 //rebinned PbPb
	static const double CentBins[NC+1] = {0,5,10,20,30,40,50,60}; //rebinned PbPb
	static const double CentBinsNat[NCNAT+1] = {0,1,2,5,10,20,30,40,50,60}; //native PbPb
	static const uint CentRebin[NC+1] = {0,3,4,5,6,7,8,9}; //native -> rebinned mapping
#endif

static const double CentBinsWde[NCWDE+1] = {0,20,60};
static const uint CentRebinWde[NCWDE+1] = {0,5,9}; //native -> wide bins

static void Correlation_A_B(double *pcorr, double *pcorr_err, double num, double nume, double den, double dene){
	*pcorr = num/den;
	*pcorr_err = TMath::Abs(*pcorr)*TMath::Sqrt(
		TMath::Power(nume/num,2)+
		TMath::Power(dene/den,2));
}

static void Correlation_A_BC(double *pcorr, double *pcorr_err, double num, double nume, double den1, double den1e, double den2, double den2e){
	*pcorr = num/(den1*den2);
	*pcorr_err = TMath::Abs(*pcorr)*TMath::Sqrt(
		TMath::Power(nume/num,2)+
		TMath::Power(den1e/den1,2)+
		TMath::Power(den2e/den2,2));
}

static void Correlation_A_sqrtB(double *pcorr, double *pcorr_err, double num, double nume, double den, double dene){
	*pcorr = num/TMath::Sqrt(TMath::Abs(den));
	*pcorr_err = TMath::Abs(*pcorr)*TMath::Sqrt(
		TMath::Power(nume/num,2)+
		TMath::Power(dene/(2.0*den),2));
}

static void Correlation_A_sqrtBC(double *pcorr, double *pcorr_err, double num, double nume, double den1, double den1e, double den2, double den2e){
	*pcorr = num/TMath::Sqrt(TMath::Abs(den1*den2));
	*pcorr_err = TMath::Abs(*pcorr)*TMath::Sqrt(
		TMath::Power(nume/num,2)+
		TMath::Power(den1e/(2.0*den1),2)+
		TMath::Power(den2e/(2.0*den2),2));
}

static void Correlation_AsqrtB(double *pcorr, double *pcorr_err, double a, double ae, double b, double be){
	double t = TMath::Sqrt(TMath::Abs(b));
	*pcorr = a*t;
	*pcorr_err = TMath::Sqrt(TMath::Power(t*ae,2)+TMath::Power(a*be/(2.0*t),2));
}

static void Correlation_ABminusC(double *pcorr, double *pcorr_err, double a, double ae, double b, double be, double c, double ce){
	*pcorr = a*b-c;
	*pcorr_err = TMath::Sqrt(TMath::Power(ae*b,2)+TMath::Power(be*a,2)
		+ce*ce);
}

static void Correlation_AminusBC(double *pcorr, double *pcorr_err, double a, double ae, double b, double be, double c, double ce){
	*pcorr = a-b*c;
	*pcorr_err = TMath::Sqrt(ae*ae
		+TMath::Power(be*c,2)+TMath::Power(ce*b,2));
}

static void Correlation_ABminusCD(double *pcorr, double *pcorr_err, double a, double ae, double b, double be, double c, double ce, double d, double de){
	*pcorr = a*b-c*d;
	*pcorr_err = TMath::Sqrt(TMath::Power(ae*b,2)+TMath::Power(be*a,2)+TMath::Power(ce*d,2)+TMath::Power(de*c,2));
}

static TH1D * LoadRebin(std::function<TH1D*(uint)> load, uint a, uint b, bool debug = false){
	static bool warned = false;
	TH1D *ph_rb = 0;
	for(uint ic = a; ic < b; ++ic){
		TH1D *pgr = load(ic);
		if(!pgr){
			if(!warned){
				printf("Warning: no graph for CentBin%02u!\n",ic);
				warned = true;
			}
			continue;
		}
		if(!ph_rb){
			if(debug)
				printf("base %u\n",ic);
			ph_rb = (TH1D*)pgr->Clone();
		}else{
			if(debug)
				printf("add %u\n",ic);
			ph_rb->Add(pgr);
		}
	}

	return ph_rb;
}


static std::tuple<double, double> GetTProfileMeanAndError(std::function<TProfile*(uint)> load, uint binLow, uint binHigh, uint ih) {
	double mean = 0.;
	double error = 0.;
	double weight = 0.;
	for (int ic = binLow; ic < binHigh; ++ic) {
		TProfile *prof = load(ic);
		if(!prof){
			printf("Warning: no graph for CentBin%02u!\n",ic);
			continue;
		}else{
			mean += prof->GetBinContent(ih+1) * prof->GetBinEffectiveEntries(ih+1);
			error += prof->GetBinError(ih+1) * prof->GetBinError(ih+1) * prof->GetBinEffectiveEntries(ih+1) * prof->GetBinEffectiveEntries(ih+1);
			weight += prof->GetBinEffectiveEntries(ih+1);
		}
	}
	double wMean = mean/weight;
	return {wMean, TMath::Sqrt(error)/weight};
}

bool graphValuesOutOfRange(TGraphErrors* gr, double minVal, double maxVal) {
    if (!gr) return false;
    bool outOfRange = false;
    for (int i = 0; i < gr->GetN(); ++i) {
        double y = gr->GetPointY(i);
        if (y < minVal) {
            gr->SetPointY(i, minVal);
            outOfRange = true;
        } else if (y > maxVal) {
            gr->SetPointY(i, maxVal);
            outOfRange = true;
        }
    }
    return outOfRange;
}

TGraphErrors* mergeBins(TGraphErrors* graph) {
	const int nPoints = 6;
	TGraphErrors* newGraph = new TGraphErrors(nPoints);
	std::vector<std::tuple<int, int>> indices = {{0, 1}, {2, 2}, {3, 4}, {5, 6}, {7, 8}, {9, 10}};
	std::array<double, nPoints> xValues = {2.5, 7.5, 15., 25., 35., 45.};
	int k = 0;
	for (const auto& [i, j] : indices) {
		if (i==j) {
			newGraph->SetPoint(k, graph->GetPointX(i), graph->GetPointY(i));
			newGraph->SetPointError(k, graph->GetErrorX(i), graph->GetErrorY(i));
		} else {
			newGraph->SetPoint(k, xValues[k], (graph->GetPointY(i) + graph->GetPointY(j)) / 2);
			newGraph->SetPointError(k, (graph->GetErrorX(i) + graph->GetErrorX(j)) / 2, (graph->GetErrorY(i) + graph->GetErrorY(j)) / 2);
		}
		k++;
	}
	return newGraph;
}

void jfflucresults_TProfile(
	const TString fname = "~/Lataukset/AnalysisResults4.root",
	const TString containerPath = "/JNLSCLHC15o_default/jfluc",
	const TString outfile = "resultsRoots/results.root",
	const TString graphPath = "/projappl/project_2003112/hydro/jetscape-analysis-new/src/dependencies/rho_MAP24.root",
	const TString param="000"){
	//void jfflucresults_nl(const TString fname = "~/Asiakirjat/alice/legotrain_FFluc.1/data_0711/JFFluc_legotrain-CF_PbPb-2734_20160708-LHC10h_AOD86_MgFpMgFm.root", const TString containerPath = "/JFFlucAODLHC10h_Default/test"){
	//#define USE_AUTOCORR
#ifdef USE_AUTOCORR
	printf("<!> Caution: using autocorrelated terms");
#endif
	ofstream errortxtfile;
	errortxtfile.open(Form("%s.txt", outfile.Data()));
	//TFile *pfile = TFile::Open(fname);
	TFile *pfile = new TFile(fname,"read");
	if(!pfile->IsOpen()){
		printf("Unable to open %s.\n",fname.Data());
		delete pfile;
		return;
	}

	double vn2_mean[NH][NK][NC]; //<vn^2>
	double vn2_mean_err[NH][NK][NC];
	double vn2vm2_mean[NH][NK][NHH][NK][NC]; //<vn^2 vm^2>
	double vn2vm2_mean_err[NH][NK][NHH][NK][NC];

	double cn_mixed_mean_QC[NH][NC]; //<<2>_n>
	double cn_mixed_mean_QC_err[NH][NC];

	double cn_mixed_mean_QC_eta14[NH][NC]; //<<2>_n>
	double cn_mixed_mean_QC_eta14_err[NH][NC];
	double cncm_mixed_mean_QC[NH][NHH][NC]; //<exp(in(phi1-phi3)+im(phi2-phi4))>
	double cncm_mixed_mean_QC_err[NH][NHH][NC];
	//
	for(uint ih = 2; ih < NH; ++ih)
		for(uint ik = 1; ik < NK; ++ik)
			//for(uint ic = 0; ic < NC; ++ic){
			for(uint r = 0; r < NC; ++r){

				auto [mean, error] = GetTProfileMeanAndError([&](uint ic)->TProfile *{
					return (TProfile*) pfile->Get(Form("%s/hvna/hvnaK%02uCentBin%02u",containerPath.Data(),ik,ic));}, CentRebin[r], CentRebin[r+1], ih);
				vn2_mean[ih][ik][r] = mean;
				printf("cent: %u, v_%u: %.4f, error: %.4f\n",r, ih, mean, error);
				vn2_mean_err[ih][ik][r] = error;

			}

	for(uint ih = 2; ih < NH; ++ih)
		for(uint ik = 1; ik < NK; ++ik)
			for(uint ihh = 2; ihh < NHH; ++ihh)
				for(uint ikk = 1; ikk < NK; ++ikk)
					//for(uint ic = 0; ic < NC; ++ic){
					for(uint r = 0; r < NC; ++r){
						TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
							return (TH1D*)pfile->Get(Form("%s/hvn_vn/hvn_vnNH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						},CentRebin[r],CentRebin[r+1]);
						//TH1D *pgr = (TH1D*)pfile->Get(Form("%s/hvn_vn/hvn_vnNH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						if(!pgr){
							vn2vm2_mean[ih][ik][ihh][ikk][r] = 0.01;
							vn2vm2_mean_err[ih][ik][ihh][ikk][r] = 0.1;
							continue;
						}
						vn2vm2_mean[ih][ik][ihh][ikk][r] = pgr->GetMean();
						vn2vm2_mean_err[ih][ik][ihh][ikk][r] = pgr->GetMeanError();

						/*pgr = LoadRebin([&](uint ic)->TH1D *{
							return (TH1D*)pfile->Get(Form("%s/hcn_cn_2c/hcn_cn_2cNH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						},CentRebin[r],CentRebin[r+1]);
						//pgr = (TH1D*)pfile->Get(Form("%s/hcn_cn_2c/hcn_cn_2cNH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						cncm_mean_QC[ih][ik][ihh][ikk][r] = pgr->GetMean();
						cncm_mean_QC_err[ih][ik][ihh][ikk][r] = pgr->GetMeanError();
						
						pgr = LoadRebin([&](uint ic)->TH1D *{
							return (TH1D*)pfile->Get(Form("%s/hcn_cn_2c_eta10/hcn_cn_2c_eta10NH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						},CentRebin[r],CentRebin[r+1]);
						//pgr = (TH1D*)pfile->Get(Form("%s/hcn_cn_2c_eta10/hcn_cn_2c_eta10NH%02uK%02uNHH%02uKK%02uCentBin%02u",containerPath.Data(),ih,ik,ihh,ikk,ic));
						cncm_mean_QCeta10[ih][ik][ihh][ikk][r] = pgr->GetMean();
						cncm_mean_QCeta10_err[ih][ik][ihh][ikk][r] = pgr->GetMeanError();*/
					}

	for(uint ih = 2; ih < NH; ++ih)
		//for(uint ic = 0; ic < NC; ++ic){
		for(uint r = 0; r < NC; ++r){
			TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form("%s/hQC_SC2p_eta10/hQC_SC2p_eta10NH%02uCentBin%02u",containerPath.Data(),ih,ic));
			},CentRebin[r],CentRebin[r+1]);
			//TH1D *pgr = (TH1D*)pfile->Get(Form("%s/hQC_SC2p/hQC_SC2pNH%02uCentBin%02u",containerPath.Data(),ih,ic));
			//TH1D *pgr = (TH1D*)pfile->Get(Form("/%s/hQC_SC2p_eta10/hQC_SC2p_eta10NH%02uCentBin%02u",containerPath.Data(),ih,ic));
			if(!pgr){
				cn_mixed_mean_QC[ih][r] = 0.01;
				cn_mixed_mean_QC_err[ih][r] = 0.1;
				continue;
			}
			cn_mixed_mean_QC[ih][r] = pgr->GetMean();
			cn_mixed_mean_QC_err[ih][r] = pgr->GetMeanError();


			auto [mean, error] = GetTProfileMeanAndError([&](uint ic)->TProfile *{
				return (TProfile*) pfile->Get(Form("%s/hQC_SC2p_eta14/hQC_SC2p_eta14CentBin%02u",containerPath.Data(),ic));}, CentRebin[r], CentRebin[r+1], ih);

			cn_mixed_mean_QC_eta14[ih][r] = mean;
			cn_mixed_mean_QC_eta14_err[ih][r] = error;

		}

	for(uint ih = 2; ih < NH; ++ih)
		// for(uint ihh = 2, ihhm = ih < NHH?ih:NHH; ihh < ihhm; ++ihh) {
		for(uint ihh = 2; ihh < 6; ++ihh) {
			//for(uint ic = 0; ic < NC; ++ic){
			printf("ih: %d, ihh: %d\n", ih, ihh);

			for(uint r = 0; r < NC; ++r){
				TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
					return (TH1D*)pfile->Get(Form("%s/hQC_SC4p/hQC_SC4pNH%02uNHH%02uCentBin%02u",containerPath.Data(),ih,ihh,ic));
				},CentRebin[r],CentRebin[r+1]);
				//TH1D *pgr = (TH1D*)pfile->Get(Form("%s/hQC_SC4p/hQC_SC4pNH%02uNHH%02uCentBin%02u",containerPath.Data(),ih,ihh,ic));
				if(!pgr){
					cncm_mixed_mean_QC[ih][ihh][r] = 0.01;
					cncm_mixed_mean_QC_err[ih][ihh][r] = 0.1;
					continue;
				}
				cncm_mixed_mean_QC[ih][ihh][r] = pgr->GetMean();
				cncm_mixed_mean_QC_err[ih][ihh][r] = pgr->GetMeanError();
			}
		}

	TGraphErrors *pgr_vn[NH], *pgr_vn_QC[NH], *pgr_vn2[NH], *pgr_vn2_QC[NH], *pgr_vn2_QC_eta14[NH], *pgr_vn_QC4[NH], *pgr_vn_QC_eta14[NH];
	for(uint ih = 0; ih < NH; ++ih){
		pgr_vn[ih] = new TGraphErrors(NC);
		pgr_vn2[ih] = new TGraphErrors(NC);
		pgr_vn_QC[ih] = new TGraphErrors(NC);
		pgr_vn2_QC[ih] = new TGraphErrors(NC);
		pgr_vn2_QC_eta14[ih] = new TGraphErrors(NC);
		pgr_vn_QC_eta14[ih] = new TGraphErrors(NC);
		pgr_vn_QC4[ih] = new TGraphErrors(NC);
		for(uint ic = 0; ic < NC; ++ic){

			pgr_vn2[ih]->SetPoint(ic, 0.5*(CentBins[ic]+CentBins[ic+1]), vn2_mean[ih][1][ic]);
			pgr_vn2[ih]->SetPointError(ic, 0., vn2_mean_err[ih][1][ic]);

			double vn;
			if (vn2_mean[ih][1][ic]<0.0) {
				errortxtfile << "v" << ih << " PC, ic" << ic << " val: " << vn2_mean[ih][1][ic] << "\n";
				pgr_vn[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]), -0.1);
				pgr_vn[ih]->SetPointError(ic,0, 1.);
			} else {
				vn = TMath::Sqrt(TMath::Abs(vn2_mean[ih][1][ic]));
				pgr_vn[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]),vn);
				pgr_vn[ih]->SetPointError(ic,0,0.5*vn2_mean_err[ih][1][ic]/vn);
			}

			pgr_vn2_QC[ih]->SetPoint(ic, 0.5*(CentBins[ic]+CentBins[ic+1]), cn_mixed_mean_QC[ih][ic]);
			pgr_vn2_QC[ih]->SetPointError(ic, 0., cn_mixed_mean_QC_err[ih][ic]);

			pgr_vn2_QC_eta14[ih]->SetPoint(ic, 0.5*(CentBins[ic]+CentBins[ic+1]), cn_mixed_mean_QC_eta14[ih][ic]);
			pgr_vn2_QC_eta14[ih]->SetPointError(ic, 0., cn_mixed_mean_QC_eta14_err[ih][ic]);

			if (cn_mixed_mean_QC[ih][ic]<0.0) {
				errortxtfile << "v" << ih << " QC, ic" << ic << " val: " << cn_mixed_mean_QC[ih][ic] << "\n";
				pgr_vn_QC[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]), -0.1);
        	                pgr_vn_QC[ih]->SetPointError(ic,0, 1.);
			} else {
				vn = TMath::Sqrt(TMath::Abs(cn_mixed_mean_QC[ih][ic]));
				pgr_vn_QC[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]),vn);
				pgr_vn_QC[ih]->SetPointError(ic,0,0.5*cn_mixed_mean_QC_err[ih][ic]/vn);
			}
			if (cn_mixed_mean_QC_eta14[ih][ic]<0.0) {
				errortxtfile << "v" << ih << " QC, ic" << ic << " val: " << cn_mixed_mean_QC_eta14[ih][ic] << "\n";
				pgr_vn_QC_eta14[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]), -0.1);
        	                pgr_vn_QC_eta14[ih]->SetPointError(ic,0, 1.);
			} else {
				vn = TMath::Sqrt(TMath::Abs(cn_mixed_mean_QC_eta14[ih][ic]));
				pgr_vn_QC_eta14[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]),vn);
				pgr_vn_QC_eta14[ih]->SetPointError(ic,0,0.5*cn_mixed_mean_QC_eta14_err[ih][ic]/vn);
			}
			double c4 = cncm_mixed_mean_QC[ih][ih][ic]-2.0*cn_mixed_mean_QC[ih][ic]*cn_mixed_mean_QC[ih][ic];
			vn = TMath::Power(TMath::Abs(-c4),0.25);
			pgr_vn_QC4[ih]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]),vn);
			pgr_vn_QC4[ih]->SetPointError(ic,0,TMath::Sqrt((cn_mixed_mean_QC[ih][ic]*cn_mixed_mean_QC[ih][ic]*cn_mixed_mean_QC_err[ih][ic]*cn_mixed_mean_QC_err[ih][ic]+cncm_mixed_mean_QC_err[ih][ih][ic]*cncm_mixed_mean_QC_err[ih][ih][ic])/TMath::Power(-c4,1.5)));
		}
	}

	//natively binned vns
	TGraphErrors *pgr_vn_nat[NH];
	TGraphErrors *pgr_vn_QC_nat[NH];
	for(uint ih = 2; ih < NH; ++ih){
		pgr_vn_nat[ih] = new TGraphErrors(NCNAT);
		pgr_vn_QC_nat[ih] = new TGraphErrors(NCNAT);
		for(uint ic = 0; ic < NCNAT; ++ic){
			TH1D *pgr = (TH1D*)pfile->Get(Form("%s/hvna/hvnaNH%02dK%02dCentBin%02d",containerPath.Data(),ih,1,ic));
			if(!pgr){
				pgr_vn_nat[ih]->SetPoint(ic,0.5*(CentBinsNat[ic]+CentBinsNat[ic+1]),0);
				pgr_vn_nat[ih]->SetPointError(ic,0,0.5);
			}else{
				double vn2_mean = pgr->GetMean();
				double vn2_mean_err = pgr->GetMeanError();
				double vn = TMath::Sqrt(TMath::Abs(vn2_mean));
				pgr_vn_nat[ih]->SetPoint(ic,0.5*(CentBinsNat[ic]+CentBinsNat[ic+1]),vn);
				pgr_vn_nat[ih]->SetPointError(ic,0,0.5*vn2_mean_err/vn);
			}

			//
			pgr = (TH1D*)pfile->Get(Form("%s/hQC_SC2p/hQC_SC2pNH%02dCentBin%02d",containerPath.Data(),ih,ic));
			if(!pgr){
				pgr_vn_QC_nat[ih]->SetPoint(ic,0.5*(CentBinsNat[ic]+CentBinsNat[ic+1]),0);
				pgr_vn_QC_nat[ih]->SetPointError(ic,0,0.5);
			}else{
				double vn2_mean = pgr->GetMean();
				double vn2_mean_err = pgr->GetMeanError();
				double vn = TMath::Sqrt(TMath::Abs(vn2_mean));
				pgr_vn_QC_nat[ih]->SetPoint(ic,0.5*(CentBinsNat[ic]+CentBinsNat[ic+1]),vn);
				pgr_vn_QC_nat[ih]->SetPointError(ic,0,0.5*vn2_mean_err/vn);
			}
		}
	}
	errortxtfile.close();
	TGraphErrors *pgr_vn_wde[NH];
	for(uint ih = 2; ih < NH; ++ih){
		pgr_vn_wde[ih] = new TGraphErrors(NCWDE);
		for(uint r = 0; r < NCWDE; ++r){
			TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form("%s/hvna/hvnaNH%02uK%02uCentBin%02u",containerPath.Data(),ih,1,ic));
			},CentRebinWde[r],CentRebinWde[r+1],false);
			if(!pgr){
				pgr_vn_wde[ih]->SetPoint(r,0.5*(CentBinsWde[r]+CentBinsWde[r+1]),0);
				pgr_vn_wde[ih]->SetPointError(r,0,0.5);
				continue;
			}
			double vn2_mean = pgr->GetMean();
			double vn2_mean_err = pgr->GetMeanError();
			double vn = TMath::Sqrt(TMath::Abs(vn2_mean));
			pgr_vn_wde[ih]->SetPoint(r,0.5*(CentBinsWde[r]+CentBinsWde[r+1]),vn);
			pgr_vn_wde[ih]->SetPointError(r,0,0.5*vn2_mean_err/vn);
		}
	}


#ifdef USE_AUTOCORR
#define CorrnV4V2s2 CorrV4V2s2
#define CorrnV5V2sV3s CorrV5V2sV3s
#define CorrnV6V3s2 CorrV6V3s2
#define CorrnV6V2s3 CorrV6V2s3
#define CorrnV7V2s2V3s CorrV7V2s2V3s
#define CorrnV8V2sV3s2 CorrV8V2sV3s2
#define CorrnV6V2sV4s CorrV6V2sV4s
#define CorrnV7V2sV5s CorrV7V2sV5s
#define CorrnV7V3sV4s CorrV7V3sV4s
#else
#define CorrnV4V2s2 CorraV4V2s2
#define CorrnV5V2sV3s CorraV5V2sV3s
#define CorrnV6V3s2 CorraV6V3s2
#define CorrnV6V2s3 CorraV6V2s3
#define CorrnV7V2s2V3s CorraV7V2s2V3s
#define CorrnV8V2sV3s2 CorraV8V2sV3s2
#define CorrnV6V2sV4s CorraV6V2sV4s
#define CorrnV7V2sV5s CorraV7V2sV5s
#define CorrnV7V3sV4s CorraV7V3sV4s
#endif

	enum Corr{
		CorrV4V2s2v22, //<V4 V2^*2 v2^2>
		CorrV4V2s2v24, //<V4 V2^*2 v2^4>
		CorrV4V2s2, //<V4 V2^*2>
		CorrV5V2sV3sv22,
		CorrV5V2sV3s,
		CorrV5V2sV3sv32,
		CorrV6V2s3,
		CorrV6V3s2,
		CorrV7V2s2V3s,
		CorraV4V2s2, //autocorrelation removed
		CorraV5V2sV3s,
		CorraV6V3s2,
		CorrNV4V4V2V2,
		CorrNV3V3V2V2,
		CorrNV5V5V2V2,
		CorrNV5V5V3V3,
		CorrNV4V4V3V3,
		CorrV8V2sV3s2,
		CorrV8V2s4,
		CorraV6V2s3,
		CorraV7V2s2V3s,
		CorraV8V2sV3s2,
		CorrV6V2sV4s,
		CorrV7V2sV5s,
		CorrV7V3sV4s,
		CorraV6V2sV4s,
		CorraV7V2sV5s,
		CorraV7V3sV4s,
		nCorr
	};

	double correlator[nCorr][NC];
	double correlator_err[nCorr][NC];

	for(uint i = 0; i < nCorr; ++i){
		//for(uint ic = 0; ic < NC; ++ic){
		for(uint r = 0; r < NC; ++r){
			auto [mean, error] = GetTProfileMeanAndError([&](uint ic)->TProfile *{
					return (TProfile*) pfile->Get(Form("%s/h_corr/h_corrCentBin%02u",containerPath.Data(),ic));}, CentRebin[r], CentRebin[r+1], i);
			correlator[i][r] = mean;
			correlator_err[i][r] = error;
		}
	}

	TGraphErrors *pgr_corr[Obs::nObs];
	for(uint i = 0; i < Obs::nObs; ++i)
		pgr_corr[i] = new TGraphErrors(NC);

	double corr, corr_err;
	for(uint ic = 0; ic < NC; ++ic){
#define point(x)\
		pgr_corr[x]->SetPoint(ic,0.5*(CentBins[ic]+CentBins[ic+1]),corr);\
		pgr_corr[x]->SetPointError(ic,0,corr_err);

		Correlation_A_BC(&corr,&corr_err,correlator[CorrV4V2s2v22][ic],correlator_err[CorrV4V2s2v22][ic],
			correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Corr0);

		Correlation_A_BC(&corr,&corr_err,vn2_mean[2][3][ic],vn2_mean_err[2][3][ic],
			vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Corr1);

		Correlation_A_BC(&corr,&corr_err,correlator[CorrV5V2sV3sv22][ic],correlator_err[CorrV5V2sV3sv22][ic],
			correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic],vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Corr2);

		//Correlation_A_BC(&corr,&corr_err,vn2vm2_mean[2][2][3][1][ic],vn2vm2_mean_err[2][2][3][1][ic],
		//	vn2vm2_mean[2][1][3][1][ic],vn2vm2_mean_err[2][1][3][1][ic],vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		Correlation_A_BC(&corr,&corr_err,vn2vm2_mean[2][2][3][1][ic],vn2vm2_mean_err[2][2][3][1][ic],
			correlator[CorrNV3V3V2V2][ic],correlator_err[CorrNV3V3V2V2][ic],vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Corr3);

		//rho_422
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],
			vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],vn2_mean[4][1][ic],vn2_mean_err[4][1][ic]);
		point(Obs::Rho42);

		//rho_523
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic],
			vn2vm2_mean[2][1][3][1][ic],vn2vm2_mean_err[2][1][3][1][ic],vn2_mean[5][1][ic],vn2_mean_err[5][1][ic]);
		point(Obs::Rho523);

		//rho_6222
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV6V2s3][ic],correlator_err[CorrnV6V2s3][ic],
			vn2_mean[2][3][ic],vn2_mean_err[2][3][ic],vn2_mean[6][1][ic],vn2_mean_err[6][1][ic]);
		point(Obs::Rho62);

		//rho_633
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV6V3s2][ic],correlator_err[CorrnV6V3s2][ic],
			vn2_mean[3][2][ic],vn2_mean_err[3][2][ic],vn2_mean[6][1][ic],vn2_mean_err[6][1][ic]);
		point(Obs::Rho63);

		//rho_624
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV6V2sV4s][ic],correlator_err[CorrnV6V2sV4s][ic],
			correlator[CorrNV4V4V2V2][ic],correlator_err[CorrNV4V4V2V2][ic],
			vn2_mean[6][1][ic],vn2_mean_err[6][1][ic]);
		point(Obs::Rho624);

		//rho_7223
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV7V2s2V3s][ic],correlator_err[CorrnV7V2s2V3s][ic],
			vn2vm2_mean[2][2][3][1][ic],vn2vm2_mean_err[2][2][3][1][ic],vn2_mean[7][1][ic],vn2_mean_err[7][1][ic]);
		point(Obs::Rho723);

		//rho_725
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV7V2sV5s][ic],correlator_err[CorrnV7V2sV5s][ic],
			vn2vm2_mean[2][1][5][1][ic],vn2vm2_mean_err[2][1][5][1][ic],
			vn2_mean[7][1][ic],vn2_mean_err[7][1][ic]);
		point(Obs::Rho725);

		//rho_734
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV7V3sV4s][ic],correlator_err[CorrnV7V3sV4s][ic],
			vn2vm2_mean[3][1][4][1][ic],vn2vm2_mean_err[3][1][4][1][ic],
			vn2_mean[7][1][ic],vn2_mean_err[7][1][ic]);
		point(Obs::Rho734);
		
		//rho_82222
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrV8V2s4][ic],correlator_err[CorrV8V2s4][ic],
			vn2_mean[2][4][ic],vn2_mean_err[2][4][ic],vn2_mean[8][1][ic],vn2_mean_err[8][1][ic]);
		point(Obs::Rho82);

		//rho_8233
		Correlation_A_sqrtBC(&corr,&corr_err,correlator[CorrnV8V2sV3s2][ic],correlator_err[CorrnV8V2sV3s2][ic],
			vn2vm2_mean[2][1][3][2][ic],vn2vm2_mean_err[2][1][3][2][ic],vn2_mean[8][1][ic],vn2_mean_err[8][1][ic]);
		point(Obs::Rho823);
		
		//chi_422
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],
			//cn_mean_QCeta10[2][2][ic],cn_mean_QCeta10_err[2][2][ic]);
			vn2_mean[2][2][ic],vn2_mean_err[2][2][ic]);
		point(Obs::Chi42);

		//chi_422QC
		//Correlation_A_B(&corr,&corr_err,correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],
		//	//cn_mean_QCeta10[2][2][ic],cn_mean_QCeta10_err[2][2][ic]);
		//	cncm_mixed_mean_QC[2][2][ic],cncm_mixed_mean_QC_err[2][2][ic]);
		//point(Obs::Chi42QC);

		//chi_523
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic],
			//cncm_mean_QCeta10[2][1][3][1][ic],cncm_mean_QCeta10_err[2][1][3][1][ic]);
			//vn2vm2_mean[2][1][3][1][ic],vn2vm2_mean_err[2][1][3][1][ic]);
			correlator[CorrNV3V3V2V2][ic],correlator_err[CorrNV3V3V2V2][ic]);
		point(Obs::Chi523);

		//chi_6222
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV6V2s3][ic],correlator_err[CorrnV6V2s3][ic],
			//cn_mean_QCeta10[2][3][ic],cn_mean_QCeta10_err[2][3][ic]);
			vn2_mean[2][3][ic],vn2_mean_err[2][3][ic]);
		point(Obs::Chi62);

		//chi_633
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV6V3s2][ic],correlator_err[CorrnV6V3s2][ic],
			//cn_mean_QCeta10[3][2][ic],cn_mean_QCeta10_err[3][2][ic]);
			//TODO: cn2_mean[3][1]
			vn2_mean[3][2][ic],vn2_mean_err[3][2][ic]);
		point(Obs::Chi63);

		//chi_624
		double tcorra, tcorra_err, tcorrb, tcorrb_err;
		Correlation_ABminusCD(&tcorra,&tcorra_err,correlator[CorrnV6V2sV4s][ic],correlator_err[CorrnV6V2sV4s][ic],vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],correlator[CorrnV6V2s3][ic],correlator_err[CorrnV6V2s3][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic]);
		Correlation_ABminusCD(&tcorrb,&tcorrb_err,vn2_mean[4][1][ic],vn2_mean_err[4][1][ic],vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic]);
		Correlation_A_BC(&corr,&corr_err,tcorra,tcorra_err,tcorrb,tcorrb_err,vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Chi624);

		//chi_7223
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV7V2s2V3s][ic],correlator_err[CorrnV7V2s2V3s][ic],
			//cncm_mean_QCeta10[2][2][3][1][ic],cncm_mean_QCeta10_err[2][2][3][1][ic]);
			vn2vm2_mean[2][2][3][1][ic],vn2vm2_mean_err[2][2][3][1][ic]);
		point(Obs::Chi723);

		//chi_725
		Correlation_ABminusCD(&tcorra,&tcorra_err,correlator[CorrnV7V2sV5s][ic],correlator_err[CorrnV7V2sV5s][ic],correlator[CorrNV3V3V2V2][ic],correlator_err[CorrNV3V3V2V2][ic],correlator[CorrnV7V2s2V3s][ic],correlator_err[CorrnV7V2s2V3s][ic],correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic]);
		Correlation_ABminusCD(&tcorrb,&tcorrb_err,vn2_mean[5][1][ic],vn2_mean_err[5][1][ic],correlator[CorrNV3V3V2V2][ic],correlator_err[CorrNV3V3V2V2][ic],correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic],correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic]);
		Correlation_A_BC(&corr,&corr_err,tcorra,tcorra_err,tcorrb,tcorrb_err,vn2_mean[2][1][ic],vn2_mean_err[2][1][ic]);
		point(Obs::Chi725);

		//chi_734
		Correlation_ABminusCD(&tcorra,&tcorra_err,correlator[CorrnV7V3sV4s][ic],correlator_err[CorrnV7V3sV4s][ic],vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],correlator[CorrnV7V2s2V3s][ic],correlator_err[CorrnV7V2s2V3s][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic]);
		Correlation_ABminusCD(&tcorrb,&tcorrb_err,vn2_mean[4][1][ic],vn2_mean_err[4][1][ic],vn2_mean[2][2][ic],vn2_mean_err[2][2][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic]);
		Correlation_A_BC(&corr,&corr_err,tcorra,tcorra_err,tcorrb,tcorrb_err,vn2_mean[3][1][ic],vn2_mean_err[3][1][ic]);
		point(Obs::Chi734);

		//chi_823
		Correlation_A_B(&corr,&corr_err,correlator[CorrnV8V2sV3s2][ic],correlator_err[CorrnV8V2sV3s2][ic],
			//cncm_mean_QCeta10[2][1][3][2][ic],cncm_mean_QCeta10_err[2][1][3][2][ic]);
			vn2vm2_mean[2][1][3][2][ic],vn2vm2_mean_err[2][1][3][2][ic]);
		point(Obs::Chi823);

		//chi_82
		Correlation_A_B(&corr,&corr_err,correlator[CorrV8V2s4][ic],correlator_err[CorrV8V2s4][ic],
			//cn_mean_QCeta10[2][4][ic],cn_mean_QCeta10_err[2][4][ic]);
			vn2_mean[2][4][ic],vn2_mean_err[2][4][ic]);
		point(Obs::Chi82);

		//v422
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV4V2s2][ic],correlator_err[CorrnV4V2s2][ic],
			vn2_mean[2][2][ic],vn2_mean_err[2][2][ic]);
		point(Obs::V42);

		//v523
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV5V2sV3s][ic],correlator_err[CorrnV5V2sV3s][ic],
			//vn2vm2_mean[2][1][3][1][ic],vn2vm2_mean_err[2][1][3][1][ic]);
			correlator[CorrNV3V3V2V2][ic],correlator_err[CorrNV3V3V2V2][ic]);
		point(Obs::V523);

		//v6222
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV6V2s3][ic],correlator_err[CorrnV6V2s3][ic],
			vn2_mean[2][3][ic],vn2_mean_err[2][3][ic]);
		//pgr_corr[Obs::Chi62]->GetPoint(ic,tcorra,tcorrb);
		//Correlation_AsqrtB(&corr,&corr_err,tcorrb,pgr_corr[Obs::Chi62]->GetErrorY(ic),vn2_mean[2][3][ic],vn2_mean_err[2][3][ic]);
		point(Obs::V62);

		//v633
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV6V3s2][ic],correlator_err[CorrnV6V3s2][ic],
			vn2_mean[3][2][ic],vn2_mean_err[3][2][ic]);
		point(Obs::V63);

		//v624
		/*double v422, chi422, chi624, chi6222;
		pgr_corr[Obs::V42]->GetPoint(ic,tcorra,v422);
		pgr_corr[Obs::Chi42]->GetPoint(ic,tcorra,chi422);
		pgr_corr[Obs::Chi62]->GetPoint(ic,tcorra,chi6222);
		pgr_corr[Obs::Chi624]->GetPoint(ic,tcorra,chi624);
		double v4l2 = vn2_mean[4][1][ic]-v422*v422;
		//double V6V2V4 = chi624*vn2_mean[2][1][ic]*v4l2+chi6222*chi422*vn2_mean[2][3][ic];*/ //identical (verified). Relation (11) in Mode coupling paper
		//note that we're plotting v624, not v62(4l)...
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV6V2sV4s][ic],correlator_err[CorrnV6V2sV4s][ic],
			correlator[CorrNV4V4V2V2][ic],correlator_err[CorrNV4V4V2V2][ic]);
		point(Obs::V624);
		//invalid!
		/*pgr_corr[Obs::Chi624]->GetPoint(ic,tcorra,tcorrb);
		Correlation_AsqrtB(&corr,&corr_err,tcorrb,pgr_corr[Obs::Chi624]->GetErrorY(ic),vn2vm2_mean[2][1][4][1][ic],vn2vm2_mean_err[2][1][4][1][ic]);
		point(Obs::V624);*/

		//v7223
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV7V2s2V3s][ic],correlator_err[CorrnV7V2s2V3s][ic],
			vn2vm2_mean[2][2][3][1][ic],vn2vm2_mean_err[2][2][3][1][ic]);
		point(Obs::V723);

		//v725
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV7V2sV5s][ic],correlator_err[CorrnV7V2sV5s][ic],
			vn2vm2_mean[2][1][5][1][ic],vn2vm2_mean_err[2][1][5][1][ic]);
		point(Obs::V725);

		//v734
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV7V3sV4s][ic],correlator_err[CorrnV7V3sV4s][ic],
			vn2vm2_mean[3][1][4][1][ic],vn2vm2_mean_err[3][1][4][1][ic]);
		point(Obs::V734);

		//v82222
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrV8V2s4][ic],correlator_err[CorrV8V2s4][ic],
			vn2_mean[2][4][ic],vn2_mean_err[2][4][ic]);
		point(Obs::V82);

		//v8233
		Correlation_A_sqrtB(&corr,&corr_err,correlator[CorrnV8V2sV3s2][ic],correlator_err[CorrnV8V2sV3s2][ic],
			vn2vm2_mean[2][1][3][2][ic],vn2vm2_mean_err[2][1][3][2][ic]);
		point(Obs::V823);

		for(uint i = 0,
			sc[12][3] = {{0,3,2},{0,4,2},{0,4,3},{0,5,2},{0,5,3},{0,6,2},{0,6,3},{0,6,4},{0,7,2},{0,7,3},{0,8,2},{0,8,3}}; i < 10; ++i){
			//SCXYYXQC
			Correlation_AminusBC(&corr,&corr_err,cncm_mixed_mean_QC[sc[i][1]][sc[i][2]][ic],cncm_mixed_mean_QC_err[sc[i][1]][sc[i][2]][ic],
				cn_mixed_mean_QC[sc[i][1]][ic],cn_mixed_mean_QC_err[sc[i][1]][ic],cn_mixed_mean_QC[sc[i][2]][ic],cn_mixed_mean_QC_err[sc[i][2]][ic]);
				//cn_mixed_mean_QC[sc[i][2]][ic],cn_mixed_mean_QC_err[sc[i][2]][ic],cn_mixed_mean_QC[sc[i][1]][ic],cn_mixed_mean_QC_err[sc[i][1]][ic]);
			point(Obs::SC3223QC+i);

			double corr1 = corr, corr_err1 = corr_err;

			//SCXYYXNQC
			Correlation_A_BC(&corr,&corr_err,corr1,corr_err1,
				vn2_mean[sc[i][2]][1][ic],vn2_mean_err[sc[i][2]][1][ic],vn2_mean[sc[i][1]][1][ic],vn2_mean_err[sc[i][1]][1][ic]);
			point(Obs::SC3223NQC+i);

			//SCXYYXNQCNG
			Correlation_A_BC(&corr,&corr_err,corr1,corr_err1,
				cn_mixed_mean_QC[sc[i][2]][ic],cn_mixed_mean_QC_err[sc[i][2]][ic],cn_mixed_mean_QC[sc[i][1]][ic],cn_mixed_mean_QC_err[sc[i][1]][ic]);
			point(Obs::SC3223NQCNG+i);
		}
	}

	const char *pslabel[] = {"Charged","Proton","Pion","Kaon"};

	TGraphErrors *pgr_ptmean[4];
	TGraphErrors *pgr_mult[4];
	TGraphErrors *pgr_mult_rap[4];
	TString spath_pT[4], spath_mult[4], spath_mult_rap[4];
	spath_pT[0].Form("%s/hChargedPtJacek/hChargedPtJacekCentBin%%02u",containerPath.Data());
	spath_mult[0].Form("%s/h_mult_eta05/h_mult_eta05CentBin%%02u",containerPath.Data());
	spath_mult_rap[0].Form("%s/h_mult_rap05/h_mult_rap05CentBin%%02u",containerPath.Data());
	for(uint is = 1; is < 4; ++is){
		spath_pT[is].Form("%s/hPtJacek_%s/hPtJacek_%sCentBin%%02u",containerPath.Data(),pslabel[is],pslabel[is]);
		spath_mult[is].Form("%s/h_mult_eta05_%s/h_mult_eta05_%sCentBin%%02u",containerPath.Data(),pslabel[is],pslabel[is]);
		spath_mult_rap[is].Form("%s/h_mult_rap05_%s/h_mult_rap05_%sCentBin%%02u",containerPath.Data(),pslabel[is],pslabel[is]);
	}
	
	for(uint is = 0; is < 4; ++is){
		pgr_ptmean[is] = new TGraphErrors(NC);
		pgr_mult[is] = new TGraphErrors(NC);
		pgr_mult_rap[is] = new TGraphErrors(NC);
		for(uint r = 0; r < NC; ++r){
			//process pT mean
			TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form(spath_pT[is].Data(),ic));
			},CentRebin[r],CentRebin[r+1]);
			double ptmean, ptmean_err;
			if(!pgr){
				ptmean = 0.0;
				ptmean_err = 1.0;
			}else{
				ptmean = pgr->GetMean();
				ptmean_err = pgr->GetMeanError();
			}
			pgr_ptmean[is]->SetPoint(r,0.5*(CentBins[r]+CentBins[r+1]),ptmean);
			pgr_ptmean[is]->SetPointError(r,0,ptmean_err);

			//process multiplicity
			pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form(spath_mult[is].Data(),ic));
			},CentRebin[r],CentRebin[r+1]);
			double mult, mult_err;
			if(!pgr){
				mult = 0.0;
				mult_err = 2e3;
			}else{
				mult = pgr->GetMean();
				mult_err = pgr->GetMeanError();
			}
			pgr_mult[is]->SetPoint(r,0.5*(CentBins[r]+CentBins[r+1]),mult);
			pgr_mult[is]->SetPointError(r,0,mult_err);


			//process multiplicity with rapidity
			pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form(spath_mult_rap[is].Data(),ic));
			},CentRebin[r],CentRebin[r+1]);
			double mult_rap, mult_rap_err;
			if(!pgr){
				mult_rap = 0.0;
				mult_rap_err = 2e3;
			}else{
				mult_rap = pgr->GetMean();
				mult_rap_err = pgr->GetMeanError();
			}
			pgr_mult_rap[is]->SetPoint(r,0.5*(CentBins[r]+CentBins[r+1]),mult_rap);
			pgr_mult_rap[is]->SetPointError(r,0,mult_rap_err);
		}
#if OO
		//pgr_mult[is] = mergeBins(pgr_mult[is]);
		//pgr_mult_rap[is] = mergeBins(pgr_mult_rap[is]);
		//if (is>0) pgr_ptmean[is] = mergeBins(pgr_ptmean[is]);
#endif
	}

	TGraphErrors *pgr_ecc[2];
	for(uint ih = 2; ih < 4; ++ih){
		pgr_ecc[ih-2] = new TGraphErrors(NC);
		for(uint r = 0; r < NC; ++r){
			TH1D *pgr = LoadRebin([&](uint ic)->TH1D *{
				return (TH1D*)pfile->Get(Form("%s/h_ecc/h_eccNH%02uCentBin%02u",containerPath.Data(),ih,ic));
			},CentRebin[r],CentRebin[r+1],false);
			if(!pgr){
				pgr_ecc[ih-2]->SetPoint(r,0.5*(CentBins[r]+CentBins[r+1]),0.0);
				pgr_ecc[ih-2]->SetPointError(r,0,0.5);
				continue;
			}
			pgr_ecc[ih-2]->SetPoint(r,0.5*(CentBins[r]+CentBins[r+1]),pgr->GetMean());
			pgr_ecc[ih-2]->SetPointError(r,0,pgr->GetMeanError());
		}
	}

	pfile->Close();
	delete pfile;

	//Writing to root file
	TFile *pfout = new TFile(outfile,"recreate");
	pfout->cd();

	//single_vn
	for(uint i = 2; i < NH; ++i){
		pgr_vn[i]->Write(Form("gr_v%u",i));
		pgr_vn2[i]->Write(Form("gr_v%u_2",i));
		pgr_vn_QC[i]->Write(Form("gr_v%u_QC",i));
		pgr_vn2_QC[i]->Write(Form("gr_v%u_2_QC",i));
		pgr_vn2_QC_eta14[i]->Write(Form("gr_v%u_2_QC_eta14",i));
		pgr_vn_QC_eta14[i]->Write(Form("gr_v%u_QC_eta14",i));
		if (i==2) pgr_vn_QC4[i] = mergeBins(pgr_vn_QC4[i]);
		pgr_vn_QC4[i]->Write(Form("gr_v%u_QC4",i));
		pgr_vn_nat[i]->Write(Form("gr_v%u_nat",i));
		pgr_vn_QC_nat[i]->Write(Form("gr_v%u_QC_nat",i));
		pgr_vn_wde[i]->Write(Form("gr_v%u_wde",i));

		delete pgr_vn[i];
		delete pgr_vn2[i];
		delete pgr_vn_QC[i];
		delete pgr_vn2_QC[i];
		delete pgr_vn2_QC_eta14[i];
		delete pgr_vn_QC4[i];
		delete pgr_vn_QC_eta14[i];
		delete pgr_vn_nat[i];
		delete pgr_vn_wde[i];
	}

	for(uint i = 0; i < Obs::nObs; ++i){
		if (i==Obs::SC3223NQCNG || i==Obs::SC3223NQC) graphValuesOutOfRange(pgr_corr[i], -.7, 1.0);
		if (i==Obs::SC4224NQCNG || i==Obs::SC4224NQC) graphValuesOutOfRange(pgr_corr[i], -.5, 3.0);
		pgr_corr[i]->Write(Obs::name[i]);
	}
	
	for(uint is = 0; is < 4; ++is){
		pgr_ptmean[is]->Write(Form("gr_pTmean_%s",pslabel[is]));
		pgr_mult[is]->Write(Form("gr_mult_%s",pslabel[is]));
		pgr_mult_rap[is]->Write(Form("gr_mult_rap_%s",pslabel[is]));
		delete pgr_ptmean[is];
		delete pgr_mult[is];
		delete pgr_mult_rap[is];
	}
	for(uint i = 2; i < 4; ++i){
		pgr_ecc[i-2]->Write(Form("gr_ecc%u",i));
		delete pgr_ecc[i-2];
	}
	
	pfout->Close();
	delete pfout;
}

