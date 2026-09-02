
namespace Obs{

enum Obs{
	Corr0,Corr1,Corr2,Corr3,
	Rho42,Rho523,Rho62,Rho63,Rho723,Rho82,Rho823,Rho624,Rho725,Rho734, //rho
	Chi42,Chi523,Chi62,Chi63,Chi723,Chi82,Chi823,Chi624,Chi725,Chi734, //chi
	V42,V523,V62,V63,V723,V82,V823,V624,V725,V734, //projected v_n{psi_m}
	SC3223,SC4224,SC4334,SC5225,SC5335, //SC with eta gap
	SC3223QC,SC4224QC,SC4334QC,SC5225QC,SC5335QC,SC6226QC,SC6336QC,SC6446QC,SC7227QC,SC7337QC,SC8228QC,SC8338QC, //SC with QC
	SC3223NQC,SC4224NQC,SC4334NQC,SC5225NQC,SC5335NQC,SC6226NQC,SC6336NQC,SC6446NQC,SC7227NQC,SC7337NQC,SC8228NQC,SC8338NQC, //normalized SC with QC
	SC3223NQCNG,SC4224NQCNG,SC4334NQCNG,SC5225NQCNG,SC5335NQCNG,SC6226NQCNG,SC6336NQCNG,SC6446NQCNG,SC7227NQCNG,SC7337NQCNG,SC8228NQCNG,SC8338NQCNG, //normalized SC with QC
	nObs
};

/*static const TString name[nObs] = {"gr_corr0","gr_corr1","gr_corr2","gr_corr3",
	"gr_rho0","gr_rho1","gr_rho2","gr_rho3","gr_rho4","gr_rho5","gr_rho6",
	"gr_chi0","gr_chi1","gr_chi2","gr_chi3","gr_chi4","gr_chi5","gr_chi6","gr_chi7","gr_chi8","gr_chi9",
	"gr_vnm0","gr_vnm1","gr_vnm2","gr_vnm3","gr_vnm4","gr_vnm5","gr_vnm6","gr_vnm7","gr_vnm8","gr_vnm9",
	"gr_sc0","gr_sc1","gr_sc2","gr_sc3","gr_sc4",
	"gr_sc0_QC","gr_sc1_QC","gr_sc2_QC","gr_sc3_QC","gr_sc4_QC",
	"gr_sc0N_QC","gr_sc1N_QC","gr_sc2N_QC","gr_sc3N_QC","gr_sc4N_QC"
};*/

static const TString name[nObs] = {"gr_corr0","gr_corr1","gr_corr2","gr_corr3",
	"gr_rho422","gr_rho523","gr_rho6222","gr_rho633","gr_rho7223","gr_rho82222","gr_rho8233","gr_rho624","gr_rho725","gr_rho734",
	"gr_chi422","gr_chi523","gr_chi6222","gr_chi633","gr_chi7223","gr_chi82222","gr_chi8233","gr_chi624","gr_chi725","gr_chi734",
	"gr_vnm422","gr_vnm523","gr_vnm6222","gr_vnm633","gr_vnm7223","gr_vnm82222","gr_vnm8233","gr_vnm624","gr_vnm725","gr_vnm734",
	"gr_sc0","gr_sc1","gr_sc2","gr_sc3","gr_sc4",
	"gr_sc32_QC","gr_sc42_QC","gr_sc43_QC","gr_sc52_QC","gr_sc53_QC","gr_sc62_QC","gr_sc63_QC","gr_sc64_QC","gr_sc72_QC","gr_sc73_QC","gr_sc82_QC","gr_sc83_QC",
	"gr_sc32N_QC","gr_sc42N_QC","gr_sc43N_QC","gr_sc52N_QC","gr_sc53N_QC","gr_sc62N_QC","gr_sc63N_QC","gr_sc64N_QC","gr_sc72N_QC","gr_sc73N_QC","gr_sc82N_QC","gr_sc83N_QC",
	"gr_sc32N_QC_nogap","gr_sc42N_QC_nogap","gr_sc43N_QC_nogap","gr_sc52N_QC_nogap","gr_sc53N_QC_nogap","gr_sc62N_QC_nogap","gr_sc63N_QC_nogap","gr_sc64N_QC_nogap","gr_sc72N_QC_nogap","gr_sc73N_QC_nogap","gr_sc82N_QC_nogap","gr_sc83N_QC_nogap",
};

}
