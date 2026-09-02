#include <vector>
#include <iostream>
#include <string>

#include "TProfile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"
TProfile* getProfile(TList* inl, const char* name);
TH1D* getHist(TList* inl, const char* name);
TH1D* getK2sub(TList* inl);
TH1D* getK3sub(TList* inl);
TH1D* getK4sub(TList* inl);
TH1D* getK2Normsub(TList* inl);
TH1D* getStdSkewsub(TList* inl);
TH1D* getIntSkewsub(TList* inl);
TH1D* getExKursub(TList* inl);
TH1D* getCk(TList* inl);
TH1D* getVarV22(TList* inl);
TH1D* getCV22pt(TList* inl);
TH1D* getRhoV22pt(TList* inl);
TH1D* getCV24pt2(TList* inl);
TH1D* getCV24pt2_term1(TList* inl);
TH1D* getCV24pt2_term2(TList* inl);
TH1D* getCV24pt2_term3(TList* inl);
TH1D* getCV24pt2_term4(TList* inl);
TH1D* getCV24pt2_term5(TList* inl);
TH1D* histPower(TH1D* inhist, double power);
template<typename T>
TGraphErrors* histToGraph(T h);
int main(int argc, const char* argv[]){
    if(argc<3) { std::cout << "No input file and/or output folder" << std::endl; return -1; }
    std::string filename = argv[1];
    std::string outputfolder = argv[2];
    TFile* f = TFile::Open(filename.c_str());
    TList* l = (TList*)f->Get("hydro/vnptCorr/OutputPtVnCorr");

    TGraphErrors* c2sub = histToGraph(getK2sub(l));
    TGraphErrors* k2sub = histToGraph(getK2Normsub(l));
    TGraphErrors* stdskewsub = histToGraph(getStdSkewsub(l));
    TGraphErrors* intskewsub = histToGraph(getIntSkewsub(l));
    TGraphErrors* exkursub = histToGraph(getExKursub(l));
    TGraphErrors* cv22pt = histToGraph(getCV22pt(l));
    TGraphErrors* ck = histToGraph(getCk(l));
    TGraphErrors* varv22 = histToGraph(getVarV22(l));
    TGraphErrors* rhov22pt = histToGraph(getRhoV22pt(l));
    TGraphErrors* cv24pt = histToGraph(getCV24pt2(l));
    TGraphErrors* cv24pt_1 = histToGraph(getCV24pt2_term1(l));
    TGraphErrors* cv24pt_2 = histToGraph(getCV24pt2_term2(l));
    TGraphErrors* cv24pt_3 = histToGraph(getCV24pt2_term3(l));
    TGraphErrors* cv24pt_4 = histToGraph(getCV24pt2_term4(l));
    TGraphErrors* cv24pt_5 = histToGraph(getCV24pt2_term5(l));

    TFile* fout = new TFile(Form("%s/output_%s", outputfolder.data(), filename.c_str()) ,"RECREATE");
    c2sub->Write();
    k2sub->Write();
    stdskewsub->Write();
    intskewsub->Write();
    exkursub->Write();
    cv22pt->Write();
    ck->Write();
    varv22->Write();
    rhov22pt->Write();
    cv24pt->Write();
/*     cv24pt_1->Write();
    cv24pt_2->Write();
    cv24pt_3->Write();
    cv24pt_4->Write();
    cv24pt_5->Write(); */
    fout->Close();

    return 0;

}
TProfile* getProfile(TList* inl, const char* name){
    if(!dynamic_cast<TProfile*>(inl->FindObject(name))) {
        std::cout << "No object of name " << name << " in list " << inl->GetName() << std::endl;
        return nullptr;
    }
    else {
        TProfile* pf = dynamic_cast<TProfile*>(inl->FindObject(name));
        pf->Sumw2();
        for(int bin = 1; bin<=pf->GetNbinsX();++bin){
            std::cout << "Setting scaled error for " << name << ": before (" << pf->GetBinError(bin) << "), after (" << pf->GetBinError(bin)/std::sqrt(pf->GetBinEntries(bin)) << ")" << std::endl;
            double be = pf->GetBinError(bin);
            double bn = std::sqrt(pf->GetBinEntries(bin));
            std::cout << "be/bn = " << be/bn << std::endl;
            pf->SetBinError(bin,be/bn);
            std::cout << "Errors set to " << pf->GetBinError(bin) << std::endl;
        }
        return pf;
    }
}
TH1D* getHist(TList* inl, const char* name){
    if(!dynamic_cast<TProfile*>(inl->FindObject(name))) {
        std::cout << "No object of name " << name << " in list " << inl->GetName() << std::endl;
	exit(1);
        return nullptr;
    }
    else {
        TProfile* pf = dynamic_cast<TProfile*>(inl->FindObject(name));
        TH1D* h = dynamic_cast<TH1D*>(pf->ProjectionX(name));
        for(int bin = 1; bin<=pf->GetNbinsX();++bin){
            std::cout << "pf error = " << pf->GetBinError(bin) << " , profile entries = " << std::sqrt(pf->GetBinEffectiveEntries(bin)) << std::endl;
            h->SetBinError(bin,pf->GetBinError(bin)/std::sqrt(pf->GetBinEffectiveEntries(bin)));
            std::cout << h->GetBinError(bin) << std::endl;
        }
        return h;
    }
}
TH1D* getK2sub(TList* inl){
    TH1D* reth = getHist(inl,"profile_pt2AB");
    TH1D* mptA = getHist(inl,"profile_ptA");
    TH1D* mptB = getHist(inl,"profile_ptB");
    mptA->Multiply(mptB);
    reth->Add(mptA,-1);
    reth->SetName("c2sub");
    return reth;
}
TH1D* getK3sub(TList* inl){
    TH1D* term1_1 = getHist(inl,"profile_pt3AAB");

    TH1D* term2_1 = getHist(inl,"profile_pt2AB");
    TH1D* mptA = getHist(inl,"profile_ptA");
    term2_1->Multiply(mptA);

    TH1D* term3_1 = getHist(inl,"profile_pt2AA");
    TH1D* mptB = getHist(inl,"profile_ptB");
    term3_1->Multiply(mptB);

    TH1D* term4_1 = dynamic_cast<TH1D*>(mptA->Clone("term4_1"));
    term4_1->Multiply(mptA);
    term4_1->Multiply(mptB);

    term1_1->Add(term2_1,-2);
    term1_1->Add(term3_1,-1);
    term1_1->Add(term4_1,2);

    //Now swap subevents to consider two particles from subevent 2
    TH1D* term1_2 = getHist(inl,"profile_pt3ABB");

    TH1D* term2_2 = getHist(inl,"profile_pt2AB");
    term2_2->Multiply(mptB);

    TH1D* term3_2 = getHist(inl,"profile_pt2BB");
    term3_2->Multiply(mptA);

    TH1D* term4_2 = dynamic_cast<TH1D*>(mptB->Clone("term4_2"));
    term4_2->Multiply(mptB);
    term4_2->Multiply(mptA);

    term1_2->Add(term2_2,-2);
    term1_2->Add(term3_2,-1);
    term1_2->Add(term4_2,2);

    //Take average of two cases
    TH1D* reth = (TH1D*)term1_2->Clone("reth_k3sub");
    reth->Add(term1_1);
    reth->Scale(0.5);
    reth->SetName("k3sub");

    return reth;
}
TH1D* getK4sub(TList* inl){
    TH1D* reth = getHist(inl,"profile_pt4AABB");

    TH1D* term2 = getHist(inl,"profile_pt3AAB");
    TH1D* mptB = getHist(inl,"profile_ptB");
    term2->Multiply(mptB);

    TH1D* term3 = getHist(inl,"profile_pt3ABB");
    TH1D* mptA = getHist(inl,"profile_ptA");
    term3->Multiply(mptA);

    TH1D* term4 = getHist(inl,"profile_pt2AA");
    TH1D* mpt2_BB = getHist(inl,"profile_pt2BB");
    term4->Multiply(mpt2_BB);

    TH1D* term5 = getHist(inl,"profile_pt2AB");
    term5->Multiply(term5);

    TH1D* term6 = getHist(inl,"profile_pt2AA");
    term6->Multiply(mptB);
    term6->Multiply(mptB);

    TH1D* term7 = getHist(inl,"profile_pt2BB");
    term7->Multiply(mptA);
    term7->Multiply(mptA);

    TH1D* term8 = getHist(inl,"profile_pt2AB");
    term8->Multiply(mptA);
    term8->Multiply(mptB);

    TH1D* term9 = dynamic_cast<TH1D*>(mptA->Clone("term9"));
    term9->Multiply(mptA);
    term9->Multiply(mptB);
    term9->Multiply(mptB);

    reth->Add(term2,-2);
    reth->Add(term3,-2);
    reth->Add(term4,-1);
    reth->Add(term5,-2);
    reth->Add(term6,2);
    reth->Add(term7,2);
    reth->Add(term8,8);
    reth->Add(term9,-6);

    reth->SetName("k4sub");
    return reth;
}
TH1D* getK2Normsub(TList* inl){
    TH1D* reth = getK2sub(inl);
    TH1D* mptA = getHist(inl,"profile_ptA");
    TH1D* mptB = getHist(inl,"profile_ptB");
    mptA->Multiply(mptB);
    reth->Divide(mptA);
    reth->SetName("k2sub");
    return reth;
}
TH1D* getStdSkewsub(TList* inl){
    TH1D* reth = getK3sub(inl);
    TH1D* k2sub = getK2sub(inl);
    reth->Divide(histPower(k2sub,3./2));
    reth->SetName("stdskewsub");
    return reth;

}
TH1D* getIntSkewsub(TList* inl){
    TH1D* reth = getK3sub(inl);
    TH1D* k2sub = getK2sub(inl);
    TH1D* mptA = getHist(inl,"profile_ptA"); //mptA
    TH1D* mptB = getHist(inl,"profile_ptB"); //mptB
    mptA->Add(mptB);
    mptA->Scale(0.5);
    reth->Multiply(mptA);
    reth->Divide(histPower(k2sub,2));
    reth->SetName("intskewsub");
    return reth;
}
TH1D* getExKursub(TList* inl){
    TH1D* reth = getK4sub(inl);
    TH1D* k2sub = getK2sub(inl);
    reth->Divide(histPower(k2sub,2));
    reth->SetName("exkursub");
    return reth;
}
TH1D* getCV22pt(TList* inl){
    TH1D* reth = getHist(inl,"profile_v22pt_gap");
    TH1D* v22 = getHist(inl,"profile_v22_gap");
    TH1D* mpt = getHist(inl,"profile_ptmid");
    mpt->Multiply(v22);
    reth->Add(mpt,-1);
    reth->SetName("cv22pt");
    return reth;
}
TH1D* getCk(TList* inl){
    TH1D* reth = getHist(inl,"profile_pt2mid");
    TH1D* mpt = getHist(inl,"profile_ptmid");
    mpt->Multiply(mpt);
    reth->Add(mpt,-1);
    reth->SetName("ck");
    return reth;
}
TH1D* getVarV22(TList* inl){
    TH1D* reth = getHist(inl,"profile_v24");
    TH1D* v22 = getHist(inl,"profile_v22");
    TH1D* v22_gap = getHist(inl,"profile_v22_gap");
    v22->Multiply(v22);
    v22_gap->Multiply(v22_gap);
    reth->Add(v22,-1);
    //reth->Add(v22_gap);
    reth->SetName("varv22");
    return reth;
}
TH1D* getRhoV22pt(TList* inl){
    TH1D* reth = getCV22pt(inl);
    TH1D* varv22 = getVarV22(inl);
    TH1D* ck = getCk(inl);
    TH1D* sqck = histPower(ck,1./2);
    TH1D* sqvarv22 = histPower(varv22,1./2);
    sqck->Multiply(sqvarv22);
    reth->Divide(sqck);
    reth->SetName("rhov22pt");
    return reth;
}
TH1D* getCV24pt2(TList* inl){
       // First get the <v2^4*deltapt^2> term
        auto reth = getHist(inl,"profile_v24pt2");
        auto v24pt_mpt = getHist(inl,"profile_v24pt_6pc_w");
        auto mpt = getHist(inl,"profile_pt");
        v24pt_mpt->Multiply(mpt);
        auto v24_mpt_mpt = getHist(inl,"profile_v24");
        v24_mpt_mpt->Multiply(mpt);
        v24_mpt_mpt->Multiply(mpt);
        reth->Add(v24pt_mpt,-2);
        reth->Add(v24_mpt_mpt);

        // Then subtract the 4 x <v2^2>*<vn^2*deltapt^2> term
        auto v22deltapt2_v22 = getHist(inl,"profile_v22pt2");
        auto v22pt_mpt = getHist(inl,"profile_v22pt_4pc_w");
        v22pt_mpt->Multiply(mpt);
        auto v22_mpt_mpt = getHist(inl,"profile_v22_4pc_w");
        v22_mpt_mpt->Multiply(mpt);
        v22_mpt_mpt->Multiply(mpt);
        v22deltapt2_v22->Add(v22pt_mpt,-2);
        v22deltapt2_v22->Add(v22_mpt_mpt);
        auto v22 = getHist(inl,"profile_v22");
        v22deltapt2_v22->Multiply(v22);
        reth->Add(v22deltapt2_v22,-4);

        // Then subtract the 4 x <v2^2*deltapt>*<v2^2*deltapt> term
        auto v22deltapt_v22deltapt = getHist(inl,"profile_v22pt");
        auto v22_mpt = getHist(inl,"profile_v22_3pc_w");
        v22_mpt->Multiply(mpt);
        v22deltapt_v22deltapt->Add(v22_mpt,-1);
        v22deltapt_v22deltapt->Multiply(v22deltapt_v22deltapt);
        reth->Add(v22deltapt_v22deltapt,-4);

        // Then add the 4 x <v2^2>*<v2^2>*<deltapt^2> term
        auto v22_v22_deltapt2 = getHist(inl,"profile_v22");
        v22_v22_deltapt2->Multiply(v22_v22_deltapt2);
        auto dpt2 = getHist(inl,"profile_pt2");
        auto mpt_mpt = dynamic_cast<TH1D*>(mpt->Clone("<mpt><mpt>"));
        mpt_mpt->Multiply(mpt);
        dpt2->Add(mpt_mpt,-1);
        v22_v22_deltapt2->Multiply(dpt2);
        reth->Add(v22_v22_deltapt2,4);

        // Then subtract the <v2^4>*<deltapt^2> term
        auto v24_deltapt2 = getHist(inl,"profile_v24");
        v24_deltapt2->Multiply(dpt2);
        reth->Add(v24_deltapt2,-1);

        reth->SetName("cv24pt2");
        return reth;
}
TH1D* getCV24pt2_term1(TList* inl){
    // First get the <v2^4*deltapt^2> term
     auto reth = getHist(inl,"profile_v24pt2");
     auto v24pt_mpt = getHist(inl,"profile_v24pt_6pc_w");
     auto mpt = getHist(inl,"profile_pt");
     v24pt_mpt->Multiply(mpt);
     auto v24_mpt_mpt = getHist(inl,"profile_v24");
     v24_mpt_mpt->Multiply(mpt);
     v24_mpt_mpt->Multiply(mpt);
     reth->Add(v24pt_mpt,-2);
     reth->Add(v24_mpt_mpt);

     reth->SetName("cv24pt2_term1");
     return reth;
}
TH1D* getCV24pt2_term2(TList* inl){
    // Then subtract the 4 x <v2^2>*<vn^2*deltapt^2> term
    auto v22deltapt2_v22 = getHist(inl,"profile_v22pt2");
    auto v22pt_mpt = getHist(inl,"profile_v22pt_4pc_w");
    auto mpt = getHist(inl,"profile_pt");
    v22pt_mpt->Multiply(mpt);
    auto v22_mpt_mpt = getHist(inl,"profile_v22_4pc_w");
    v22_mpt_mpt->Multiply(mpt);
    v22_mpt_mpt->Multiply(mpt);
    v22deltapt2_v22->Add(v22pt_mpt,-2);
    v22deltapt2_v22->Add(v22_mpt_mpt);
    auto v22 = getHist(inl,"profile_v22");
    v22deltapt2_v22->Multiply(v22);
    v22deltapt2_v22->Scale(-4);
    v22deltapt2_v22->SetName("cv24pt2_term2");
    return v22deltapt2_v22;
}
TH1D* getCV24pt2_term3(TList* inl){
    // Then subtract the 4 x <v2^2*deltapt>*<v2^2*deltapt> term
    auto v22deltapt_v22deltapt = getHist(inl,"profile_v22pt");
    auto v22_mpt = getHist(inl,"profile_v22_3pc_w");
    auto mpt = getHist(inl,"profile_pt");
    v22_mpt->Multiply(mpt);
    v22deltapt_v22deltapt->Add(v22_mpt,-1);
    v22deltapt_v22deltapt->Multiply(v22deltapt_v22deltapt);

    v22deltapt_v22deltapt->Scale(-4);
    v22deltapt_v22deltapt->SetName("cv24pt2_term3");
    return v22deltapt_v22deltapt;
}
TH1D* getCV24pt2_term4(TList* inl){
    // Then add the 4 x <v2^2>*<v2^2>*<deltapt^2> term
    auto v22_v22_deltapt2 = getHist(inl,"profile_v22");
    v22_v22_deltapt2->Multiply(v22_v22_deltapt2);
    auto mpt = getHist(inl,"profile_pt");
    auto dpt2 = getHist(inl,"profile_pt2");
    auto mpt_mpt = dynamic_cast<TH1D*>(mpt->Clone("<mpt><mpt>"));
    mpt_mpt->Multiply(mpt);
    dpt2->Add(mpt_mpt,-1);
    v22_v22_deltapt2->Multiply(dpt2);

    v22_v22_deltapt2->Scale(4);
    v22_v22_deltapt2->SetName("cv24pt2_term4");
    return v22_v22_deltapt2;
}
TH1D* getCV24pt2_term5(TList* inl){
    // Then subtract the <v2^4>*<deltapt^2> term
    auto v24_deltapt2 = getHist(inl,"profile_v24");
    auto mpt = getHist(inl,"profile_pt");
    auto dpt2 = getHist(inl,"profile_pt2");
    auto mpt_mpt = dynamic_cast<TH1D*>(mpt->Clone("<mpt><mpt>"));
    mpt_mpt->Multiply(mpt);
    dpt2->Add(mpt_mpt,-1);
    v24_deltapt2->Multiply(dpt2);

    v24_deltapt2->Scale(-1);
    v24_deltapt2->SetName("cv24pt2_term5");
    return v24_deltapt2;
}

TH1D* histPower(TH1D* inhist, double power){
    TH1D* outhist = dynamic_cast<TH1D*>(inhist->Clone(Form("power%.1f%s",power,inhist->GetName())));
    for(int bin = 1; bin<=inhist->GetNbinsX();++bin){
    if(inhist->GetBinContent(bin)>=0){
        outhist->SetBinContent(bin,std::pow(inhist->GetBinContent(bin),power));
        outhist->SetBinError(bin,power*std::pow(inhist->GetBinContent(bin),power-1)*inhist->GetBinError(bin));
    }
    else
    {
        outhist->SetBinContent(bin,-999);
        outhist->SetBinError(bin,0.000000001);
    }
    }
    return outhist;
}
template<typename T>
TGraphErrors* histToGraph(T h){
    int nBins = h->GetNbinsX();
    std::vector<double> binEdges;
    for(auto i(0); i<=nBins; ++i) binEdges.push_back(h->GetBinLowEdge(i+1));
    std::vector<double> val;
    std::vector<double> ey;
    std::vector<double> ex;
    std::vector<double> x;
    TGraphErrors* gr;
    for(Int_t bin(0); bin < nBins; ++bin)
    {
        val.push_back(h->GetBinContent(bin+1));
        ey.push_back(h->GetBinError(bin+1));
        ex.push_back(h->GetBinWidth(bin+1)*0.18);
        x.push_back((binEdges[bin]+binEdges[bin+1])/2);
    }
    gr = new TGraphErrors(nBins,&x[0],&val[0],&ex[0],&ey[0]);
    gr->SetFillStyle(1001);
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetFillColorAlpha(h->GetMarkerColor(),0.4);
    gr->SetLineColor(h->GetLineColor());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetName(h->GetName());
    return gr;
}
