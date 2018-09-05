#include <TSystem.h>
#include <TMath.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TPaveStats.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <string>
#include "SONGKYO.h"
#include "commonUtility.h"
#include "SetHLT.h"

bool usetnp = true;

void getCorrectedEffErr(const int nbins, TH1D *hrec, TH1D *hgen, TH1D *heff) {
  for (int a=0; a<nbins; a++) {
    double genInt = hgen->GetBinContent(a+1);
    double genErr = hgen->GetBinError(a+1);
    double recInt = hrec->GetBinContent(a+1);
    double recErr = hrec->GetBinError(a+1);
    double eff = recInt / genInt;

    double tmpErrGen1 = TMath::Power(eff,2) / TMath::Power(genInt,2);
    double tmpErrRec1 = TMath::Power(recErr,2);
    double tmpErr1 = tmpErrGen1 * tmpErrRec1;

    double tmpErrGen2 = TMath::Power(1-eff,2) / TMath::Power(genInt,2);
    double tmpErrRec2 = TMath::Abs(TMath::Power(genErr,2) - TMath::Power(recErr,2));
    double tmpErr2 = tmpErrGen2 * tmpErrRec2;
    double effErr = TMath::Sqrt(tmpErr1 + tmpErr2);

    if (genInt == 0) {
      heff->SetBinContent(a+1, 0);
      heff->SetBinError(a+1, 0);
    } else {
      heff->SetBinContent(a+1, eff);
      heff->SetBinError(a+1, effErr);
    }
  }
}

bool AcceptanceCut (TLorentzVector* Muon)
{
  return ((fabs(Muon->Eta())<1.2&&Muon->Pt()>=3.3)||(1.2<=fabs(Muon->Eta())&&fabs(Muon->Eta())<2.1&&Muon->Pt()>=3.9-fabs(Muon->Eta()))||(2.1<=fabs(Muon->Eta())&&Muon->Pt()>=1.3));
};

//main
void getEff_Jared_SingleMu_muGun(TString tagString = "HLT_HIL3Mu3_v1", int kTagPath=1){

  TH1::SetDefaultSumw2();

  TChain *fcha = new TChain("hionia/myTree");
  //TChain *fcha_hltBit = new TChain("hltbitanalysis/HltTree");
/*
  fcha->Add("../Trees_v1/Ups1SMM_0_30_OniaForest1.root");//PromptReco_Double_All(8TeV)
  //fcha->Add("../Trees_v1/Ups1SMM_0_30_OniaForest2.root");//PromptReco_Double_All(8TeV)
  //fcha->Add("../Trees_v1/Ups1SMM_0_30_OniaForest3.root");//PromptReco_Double_All(8TeV)
  //fcha->Add("../Trees_v1/Ups1SMM_0_30_OniaForest4.root");//PromptReco_Double_All(8TeV)
  fcha_hltBit->Add("../Trees_v1/Ups1SMM_0_30_OniaForest1.root");//PromptReco_Double_All(8TeV)
  //fcha_hltBit->Add("../Trees_v1/Ups1SMM_0_30_OniaForest2.root");//PromptReco_Double_All(8TeV)
  //fcha_hltBit->Add("../Trees_v1/Ups1SMM_0_30_OniaForest3.root");//PromptReco_Double_All(8TeV)
  //fcha_hltBit->Add("../Trees_v1/Ups1SMM_0_30_OniaForest4.root");//PromptReco_Double_All(8TeV)
*/
  fcha->Add("/eos/cms/store/group/phys_heavyions/dileptons/pshukla/OniaTrigger18/SingleMuPt_0_30_OniaForest.root");//PromptReco_Double_All(8TeV)
  //fcha_hltBit->Add("/eos/cms/store/group/phys_heavyions/dileptons/pshukla/OniaTrigger18/SingleMuPt_0_30_OniaForest.root");//PromptReco_Double_All(8TeV)

  double ptmin = 0;
  double ptmax = 15;
  double etamin = -2.4;
  double etamax = 2.4;
  double massmin = 0;
  double massmax = 30;
  double himassmin = 70;
  double himassmax = 110;
  string date="1124";
  string dataset="Data";
  //string dataset="MC";
  //string vername="test";
  //string vername="ExpressStream_8TeV_Double_285480_285549_MB";
  string vername="PromptReco_285480-285517_TnP_v1";

 //define Trigger
 // const int Ntrig = 12;
  const int Ntrig = trigname.size();

  for(int i=0; i < trigname.size(); i++){
    if(trigname[i] == tagString.Data()) {kTagPath=i;break;}
  }

  //kTagPath = 9;
  //k = kTagPath;
  kStart = 0;
  //Set Branch
  Int_t           Centrality;
  UInt_t          runNb;
  //Int_t           Gen_mu_size; //for MC
  Int_t           Reco_mu_size;
  Int_t           Reco_QQ_size;
  Int_t           trigPrescale[1000];
  //TClonesArray    *Gen_mu_4mom; //for MC
  TClonesArray    *Reco_mu_4mom;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_mu_trig[1000];
  ULong64_t       Reco_QQ_trig[1000];
  ULong64_t       Reco_QQ_mupl_trig[1000];
  ULong64_t       Reco_QQ_mumi_trig[1000];
  ULong64_t       Reco_QQ_sign[1000];
  ULong64_t       HLTriggers;
  Bool_t          Reco_mu_isGoodMuon[1000];
  Int_t           Reco_mu_nPixWMea[1000];
  Int_t           Reco_mu_nTrkWMea[1000];
  Float_t         Reco_mu_dxy[1000];
  Float_t         Reco_mu_dz[1000];
  ULong64_t       Reco_QQ_VtxProb[1000];
  Bool_t          Reco_QQ_mumi_isGoodMuon[1000];
  Int_t           Reco_QQ_mumi_nPixWMea[1000];
  Int_t           Reco_QQ_mumi_nTrkWMea[1000];
  Float_t         Reco_QQ_mumi_dxy[1000];
  Float_t         Reco_QQ_mumi_dz[1000];
  Bool_t          Reco_QQ_mupl_isGoodMuon[1000];
  Int_t           Reco_QQ_mupl_nPixWMea[1000];
  Int_t           Reco_QQ_mupl_nTrkWMea[1000];
  Float_t         Reco_QQ_mupl_dxy[1000];
  Float_t         Reco_QQ_mupl_dz[1000];
  Int_t           Reco_QQ_mupl_SelectionType[1000];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_SelectionType[1000];   //[Reco_QQ_size]

  TBranch        *b_Centrality;
  TBranch        *b_runNb;
  TBranch        *b_trigPrescale;
  //TBranch        *b_Gen_mu_size; //for MC
  TBranch        *b_Reco_mu_size;
  TBranch        *b_Reco_QQ_size;
  //TBranch        *b_Gen_mu_4mom; //for MC
  TBranch        *b_Reco_mu_4mom;
  TBranch        *b_Reco_QQ_4mom;
  TBranch        *b_Reco_QQ_VtxProb;
  TBranch        *b_Reco_QQ_mupl_4mom;
  TBranch        *b_Reco_QQ_mumi_4mom;
  TBranch        *b_Reco_mu_trig;
  TBranch        *b_Reco_QQ_trig;
  TBranch        *b_Reco_QQ_sign;
  TBranch        *b_HLTriggers;
  TBranch        *b_Reco_mu_isGoodMuon;
  TBranch        *b_Reco_mu_nPixWMea;
  TBranch        *b_Reco_mu_nTrkWMea;
  TBranch        *b_Reco_mu_dxy;
  TBranch        *b_Reco_mu_dz;
  TBranch        *b_Reco_QQ_mupl_trig;
  TBranch        *b_Reco_QQ_mumi_trig;
  TBranch        *b_Reco_QQ_mumi_isGoodMuon;
  TBranch        *b_Reco_QQ_mumi_nPixWMea;
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;
  TBranch        *b_Reco_QQ_mumi_dxy;
  TBranch        *b_Reco_QQ_mumi_dz;
  TBranch        *b_Reco_QQ_mupl_isGoodMuon;
  TBranch        *b_Reco_QQ_mupl_nPixWMea;
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;
  TBranch        *b_Reco_QQ_mupl_dxy;
  TBranch        *b_Reco_QQ_mupl_dz;
  TBranch        *b_Reco_QQ_mupl_SelectionType;   //!
  TBranch        *b_Reco_QQ_mumi_SelectionType;   //!

  //  Gen_mu_4mom = 0; //for MC
  Reco_mu_4mom = 0;
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

  fcha->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  fcha->SetBranchAddress("runNb", &runNb, &b_runNb);
  fcha->SetBranchAddress("trigPrescale", &trigPrescale, &b_trigPrescale);
  //  fcha->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size); //for MC
  //  fcha->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom); //for MC
  fcha->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  fcha->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  fcha->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
  fcha->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  fcha->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  fcha->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  fcha->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  fcha->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  fcha->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  fcha->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  fcha->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
  fcha->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
  fcha->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  fcha->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
  fcha->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  fcha->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  fcha->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  fcha->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  fcha->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
  fcha->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  fcha->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  fcha->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  fcha->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  fcha->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
  fcha->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  fcha->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  fcha->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  fcha->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  fcha->SetBranchAddress("Reco_QQ_mupl_SelectionType", Reco_QQ_mupl_SelectionType, &b_Reco_QQ_mupl_SelectionType);
  fcha->SetBranchAddress("Reco_QQ_mumi_SelectionType", Reco_QQ_mumi_SelectionType, &b_Reco_QQ_mumi_SelectionType);

  
  Int_t           NL1Stage2Muon;
  Float_t         L1Stage2MuonPt[200];
  Float_t         L1Stage2MuonEta[200];
  Float_t         L1Stage2MuonPhi[200];
  Int_t           L1Stage2MuonChg[200];
  Int_t           NL1Stage2DiMuon;
  Float_t         L1Stage2DiMuonPt[200];
  Float_t         L1Stage2DiMuonRap[200];
  Float_t         L1Stage2DiMuonInvMass[200];
  Float_t         L1Stage2DiMuondR[200];
  Float_t         L1Stage2DiMuonChg[200];

  Int_t           HLT_HIL1DoubleMu0_v1;
  Int_t           HLT_HIL2DoubleMu0_v1;
  Int_t           HLT_HIL3DoubleMu0_v2;
  Int_t           HLT_L1SingleMuOpen_v1;
  Int_t           HLT_L1SingleMu3_v1;
  Int_t           HLT_L1SingleMu5_v1;
  Int_t           HLT_L1SingleMu7_v1;
  Int_t           HLT_L1SingleMu12_v1;
  Int_t           HLT_L1SingleMu15_v1;
  Int_t           HLT_L1SingleMu20_v1;
  Int_t           HLT_HIL1Mu3_v1;
  Int_t           HLT_HIL2Mu3_v1;
  Int_t           HLT_HIL3Mu3_v1;
  Int_t           HLT_HIL1Mu5_v1;
  Int_t           HLT_HIL2Mu5_v1;
  Int_t           HLT_HIL3Mu5_v1;
  Int_t           HLT_HIL1Mu7_v1;
  Int_t           HLT_HIL2Mu7_v1;
  Int_t           HLT_HIL3Mu7_v1;
  Int_t           HLT_HIL1Mu12_v1;
  Int_t           HLT_HIL2Mu12_v1;
  Int_t           HLT_HIL3Mu12_v1;
  Int_t           HLT_HIL1Mu15_v1;
  Int_t           HLT_HIL2Mu15_v1;
  Int_t           HLT_HIL3Mu15_v1;
  Int_t           HLT_HIL1Mu20_v1;
  Int_t           HLT_HIL2Mu20_v1;
  Int_t           HLT_HIL3Mu20_v1;
  Int_t           HLT_HIL3Mu7_v2;
  Int_t           L1_DoubleMu0_Final;
  Int_t           L1_DoubleMu0_SQ_Final;
  Int_t           L1_DoubleMu0_SQ_OS_Final;
  Int_t           L1_DoubleMuOpen_Final;
  Int_t           L1_DoubleMu10_Final;
  Int_t           L1_SingleMuOpen_Final;
  Int_t           L1_SingleMu3_Final;
  Int_t           L1_SingleMu5_Final;
  Int_t           L1_SingleMu7_Final;

  TBranch        *b_NL1Stage2Muon;
  TBranch        *b_L1Stage2MuonPt;
  TBranch        *b_L1Stage2MuonEta;
  TBranch        *b_L1Stage2MuonPhi;
  TBranch        *b_L1Stage2MuonChg;
  TBranch        *b_NL1Stage2DiMuon;
  TBranch        *b_L1Stage2DiMuonPt;
  TBranch        *b_L1Stage2DiMuonRap;
  TBranch        *b_L1Stage2DiMuonInvMass;
  TBranch        *b_L1Stage2DiMuondR;
  TBranch        *b_L1Stage2DiMuonChg;
  TBranch        *b_HLT_HIL1DoubleMuOpen_v1;
  TBranch        *b_HLT_HIL1DoubleMu0_v1;
  TBranch        *b_HLT_HIL2DoubleMu0_v1;
  TBranch        *b_HLT_HIL3DoubleMu0_v2;
  TBranch        *b_HLT_L1SingleMuOpen_v1;
  TBranch        *b_HLT_L1SingleMu3_v1;
  TBranch        *b_HLT_L1SingleMu5_v1;
  TBranch        *b_HLT_L1SingleMu7_v1;
  TBranch        *b_HLT_L1SingleMu12_v1;
  TBranch        *b_HLT_L1SingleMu15_v1;
  TBranch        *b_HLT_L1SingleMu20_v1;
  TBranch        *b_HLT_HIL1Mu3_v1;
  TBranch        *b_HLT_HIL2Mu3_v1;
  TBranch        *b_HLT_HIL3Mu3_v1;
  TBranch        *b_HLT_HIL1Mu5_v1;
  TBranch        *b_HLT_HIL2Mu5_v1;
  TBranch        *b_HLT_HIL3Mu5_v1;
  TBranch        *b_HLT_HIL1Mu7_v1;
  TBranch        *b_HLT_HIL2Mu7_v1;
  TBranch        *b_HLT_HIL3Mu7_v1;
  TBranch        *b_HLT_HIL1Mu12_v1;
  TBranch        *b_HLT_HIL2Mu12_v1;
  TBranch        *b_HLT_HIL3Mu12_v1;
  TBranch        *b_HLT_HIL1Mu15_v1;
  TBranch        *b_HLT_HIL2Mu15_v1;
  TBranch        *b_HLT_HIL3Mu15_v1;
  TBranch        *b_HLT_HIL1Mu20_v1;
  TBranch        *b_HLT_HIL2Mu20_v1;
  TBranch        *b_HLT_HIL3Mu20_v1;
  TBranch        *b_HLT_HIL3Mu7_v2;
  TBranch        *b_L1_DoubleMuOpen_Final;
  TBranch        *b_L1_DoubleMu0_Final;
  TBranch        *b_L1_DoubleMu10_Final;
  TBranch        *b_L1_DoubleMu0_SQ_Final;
  TBranch        *b_L1_DoubleMu0_SQ_OS_Final;
  TBranch        *b_L1_SingleMuOpen_Final;
  TBranch        *b_L1_SingleMu3_Final;
  TBranch        *b_L1_SingleMu5_Final;
  TBranch        *b_L1_SingleMu7_Final;
/* 
   fcha_hltBit->SetBranchAddress("NL1Stage2Muon", &NL1Stage2Muon, &b_NL1Stage2Muon); 
   fcha_hltBit->SetBranchAddress("L1Stage2MuonPt", &L1Stage2MuonPt, &b_L1Stage2MuonPt); 
   fcha_hltBit->SetBranchAddress("L1Stage2MuonEta", &L1Stage2MuonEta, &b_L1Stage2MuonEta); 
   fcha_hltBit->SetBranchAddress("L1Stage2MuonPhi", &L1Stage2MuonPhi, &b_L1Stage2MuonPhi); 
   fcha_hltBit->SetBranchAddress("L1Stage2MuonChg", &L1Stage2MuonChg, &b_L1Stage2MuonChg); 
   fcha_hltBit->SetBranchAddress("NL1Stage2DiMuon", &NL1Stage2DiMuon, &b_NL1Stage2DiMuon); 
   fcha_hltBit->SetBranchAddress("L1Stage2DiMuonPt", &L1Stage2DiMuonPt, &b_L1Stage2DiMuonPt); 
   fcha_hltBit->SetBranchAddress("L1Stage2DiMuonRap", &L1Stage2DiMuonRap, &b_L1Stage2DiMuonRap); 
   fcha_hltBit->SetBranchAddress("L1Stage2DiMuonInvMass", &L1Stage2DiMuonInvMass, &b_L1Stage2DiMuonInvMass); 
   fcha_hltBit->SetBranchAddress("L1Stage2DiMuondR", &L1Stage2DiMuondR, &b_L1Stage2DiMuondR); 
   fcha_hltBit->SetBranchAddress("L1Stage2DiMuonChg", &L1Stage2DiMuonChg, &b_L1Stage2DiMuonChg); 
   fcha_hltBit->SetBranchAddress("HLT_HIL1DoubleMu0_v1", &HLT_HIL1DoubleMu0_v1, &b_HLT_HIL1DoubleMu0_v1); 
   fcha_hltBit->SetBranchAddress("HLT_HIL2DoubleMu0_v1", &HLT_HIL2DoubleMu0_v1, &b_HLT_HIL2DoubleMu0_v1); 
   fcha_hltBit->SetBranchAddress("HLT_HIL3DoubleMu0_v2", &HLT_HIL3DoubleMu0_v2, &b_HLT_HIL3DoubleMu0_v2); 
   fcha_hltBit->SetBranchAddress("HLT_L1SingleMuOpen_v1", &HLT_L1SingleMuOpen_v1, &b_HLT_L1SingleMuOpen_v1); 
   fcha_hltBit->SetBranchAddress("HLT_L1SingleMu3_v1", &HLT_L1SingleMu3_v1, &b_HLT_L1SingleMu3_v1); 
   fcha_hltBit->SetBranchAddress("HLT_L1SingleMu5_v1", &HLT_L1SingleMu5_v1, &b_HLT_L1SingleMu5_v1); 
   fcha_hltBit->SetBranchAddress("HLT_L1SingleMu7_v1", &HLT_L1SingleMu7_v1, &b_HLT_L1SingleMu7_v1); 
   fcha_hltBit->SetBranchAddress("HLT_HIL3Mu3_v1", &HLT_HIL3Mu3_v1, &b_HLT_HIL3Mu3_v1); 
   fcha_hltBit->SetBranchAddress("HLT_HIL3Mu5_v1", &HLT_HIL3Mu5_v1, &b_HLT_HIL3Mu5_v1); 
   fcha_hltBit->SetBranchAddress("HLT_HIL3Mu7_v2", &HLT_HIL3Mu7_v2, &b_HLT_HIL3Mu7_v2); 
   fcha_hltBit->SetBranchAddress("L1_DoubleMu0_Final", &L1_DoubleMu0_Final, &b_L1_DoubleMu0_Final); 
   fcha_hltBit->SetBranchAddress("L1_DoubleMu0_SQ_Final", &L1_DoubleMu0_SQ_Final, &b_L1_DoubleMu0_SQ_Final); 
   fcha_hltBit->SetBranchAddress("L1_DoubleMu0_SQ_OS_Final", &L1_DoubleMu0_SQ_OS_Final, &b_L1_DoubleMu0_SQ_OS_Final); 
   fcha_hltBit->SetBranchAddress("L1_DoubleMuOpen_Final", &L1_DoubleMuOpen_Final, &b_L1_DoubleMuOpen_Final); 
   fcha_hltBit->SetBranchAddress("L1_DoubleMu10_Final", &L1_DoubleMu10_Final, &b_L1_DoubleMu10_Final); 
   fcha_hltBit->SetBranchAddress("L1_SingleMuOpen_Final", &L1_SingleMuOpen_Final, &b_L1_SingleMuOpen_Final); 
   fcha_hltBit->SetBranchAddress("L1_SingleMu3_Final", &L1_SingleMu3_Final, &b_L1_SingleMu3_Final); 
   fcha_hltBit->SetBranchAddress("L1_SingleMu5_Final", &L1_SingleMu5_Final, &b_L1_SingleMu5_Final); 
   fcha_hltBit->SetBranchAddress("L1_SingleMu7_Final", &L1_SingleMu7_Final, &b_L1_SingleMu7_Final); 
*/

  double ptarr[] = {0, 1,2,3, 4, 5, 6, 7, 8,9,10,11,12,13,14,15 };
  //double ptarr[] = {0, 1,2,3, 5, 7, 8, 9, 10,12,15,18,25, 40};
  //double ptarr[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,30};
  const int Nptarr = sizeof(ptarr)/sizeof(double);

  double raparr[] = {-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4};
  const int Nraparr = sizeof(raparr)/sizeof(double);

  double phiarr[] = {-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,
    0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1};
  const int Nphiarr = sizeof(phiarr)/sizeof(double);

  double centarr[] = {0, 20, 40, 60, 80, 100};
  const int Ncentarr = sizeof(centarr)/sizeof(double);

  double legmin[]={0.950,0.925,0.900,0.875,0.850,0.825,0.800,0.775,0.750,0.725,0.700};
  const int Nlegmin = sizeof(legmin)/sizeof(double);
  double legmax[]={0.925,0.900,0.875,0.850,0.825,0.800,0.775,0.750,0.725,0.700,0.675};
  const int Nlegmax = sizeof(legmax)/sizeof(double);

  double l1_dR_cut = 0.3;
  double l1_mass_cut = 3;


  //Define histograms
  TH1D *test_p[Ntrig];
  TH1D *nume_p[Ntrig];
  TH1D *nume_e[Ntrig];
  TH1D *nume_phi[Ntrig];
  TH1D *deno_p[Ntrig];
  TH1D *deno_e[Ntrig];
/*  TH1D *deno_phi[Ntrig];
  
  TH1D *nume_p_l1SQ[Ntrig];
  TH1D *nume_e_l1SQ[Ntrig];
  TH1D *deno_p_l1SQ[Ntrig];
  TH1D *deno_e_l1SQ[Ntrig];

  TH1D *nume_p_l1dR[Ntrig];
  TH1D *nume_e_l1dR[Ntrig];
  TH1D *deno_p_l1dR[Ntrig];
  TH1D *deno_e_l1dR[Ntrig];
  
  TH1D *nume_p_l1OS[Ntrig];
  TH1D *nume_e_l1OS[Ntrig];
  TH1D *deno_p_l1OS[Ntrig];
  TH1D *deno_e_l1OS[Ntrig];
  
  TH1D *nume_p_l1open[Ntrig];
  TH1D *nume_e_l1open[Ntrig];
  TH1D *deno_p_l1open[Ntrig];
  TH1D *deno_e_l1open[Ntrig];
  
  TH1D *nume_p_l1SQOS[Ntrig];
  TH1D *nume_e_l1SQOS[Ntrig];
  TH1D *deno_p_l1SQOS[Ntrig];
  TH1D *deno_e_l1SQOS[Ntrig];
  
  TH1D *nume_p_l1mass[Ntrig];
  TH1D *nume_e_l1mass[Ntrig];
  TH1D *deno_p_l1mass[Ntrig];
  TH1D *deno_e_l1mass[Ntrig];
*/  
  TGraphAsymmErrors *eff_p[Ntrig];
  TGraphAsymmErrors *eff_e[Ntrig];
/*  TGraphAsymmErrors *eff_p_l1SQ[Ntrig];
  TGraphAsymmErrors *eff_e_l1SQ[Ntrig];
  TGraphAsymmErrors *eff_p_l1OS[Ntrig];
  TGraphAsymmErrors *eff_e_l1OS[Ntrig];
  TGraphAsymmErrors *eff_p_l1SQOS[Ntrig];
  TGraphAsymmErrors *eff_e_l1SQOS[Ntrig];
  TGraphAsymmErrors *eff_p_l1dR[Ntrig];
  TGraphAsymmErrors *eff_e_l1dR[Ntrig];
  TGraphAsymmErrors *eff_p_l1mass[Ntrig];
  TGraphAsymmErrors *eff_e_l1mass[Ntrig];
  TGraphAsymmErrors *eff_p_l1open[Ntrig];
  TGraphAsymmErrors *eff_e_l1open[Ntrig];
*/
  //make histrograms
  for(int i=0; i<Ntrig; i++){
    nume_p[i] = new TH1D(Form("nume_p%d",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e[i] = new TH1D(Form("nume_e%d",i),";#eta;Events",Nraparr-1,raparr);
    deno_p[i] = new TH1D(Form("deno_p%d",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e[i] = new TH1D(Form("deno_e%d",i),";#eta;Events",Nraparr-1,raparr);
/*
    nume_p_l1SQ[i] = new TH1D(Form("nume_p%d_l1SQ",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1SQ[i] = new TH1D(Form("nume_e%d_l1SQ",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1SQ[i] = new TH1D(Form("deno_p%d_l1SQ",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1SQ[i] = new TH1D(Form("deno_e%d_l1SQ",i),";#eta;Events",Nraparr-1,raparr);

    nume_p_l1dR[i] = new TH1D(Form("nume_p%d_l1dR",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1dR[i] = new TH1D(Form("nume_e%d_l1dR",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1dR[i] = new TH1D(Form("deno_p%d_l1dR",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1dR[i] = new TH1D(Form("deno_e%d_l1dR",i),";#eta;Events",Nraparr-1,raparr);

    nume_p_l1SQOS[i] = new TH1D(Form("nume_p%d_l1SQOS",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1SQOS[i] = new TH1D(Form("nume_e%d_l1SQOS",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1SQOS[i] = new TH1D(Form("deno_p%d_l1SQOS",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1SQOS[i] = new TH1D(Form("deno_e%d_l1SQOS",i),";#eta;Events",Nraparr-1,raparr);

    nume_p_l1mass[i] = new TH1D(Form("nume_p%d_l1mass",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1mass[i] = new TH1D(Form("nume_e%d_l1mass",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1mass[i] = new TH1D(Form("deno_p%d_l1mass",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1mass[i] = new TH1D(Form("deno_e%d_l1mass",i),";#eta;Events",Nraparr-1,raparr);

    nume_p_l1open[i] = new TH1D(Form("nume_p%d_l1open",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1open[i] = new TH1D(Form("nume_e%d_l1open",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1open[i] = new TH1D(Form("deno_p%d_l1open",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1open[i] = new TH1D(Form("deno_e%d_l1open",i),";#eta;Events",Nraparr-1,raparr);

    nume_p_l1OS[i] = new TH1D(Form("nume_p%d_l1OS",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e_l1OS[i] = new TH1D(Form("nume_e%d_l1OS",i),";#eta;Events",Nraparr-1,raparr);
    deno_p_l1OS[i] = new TH1D(Form("deno_p%d_l1OS",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e_l1OS[i] = new TH1D(Form("deno_e%d_l1OS",i),";#eta;Events",Nraparr-1,raparr);

    nume_phi[i] = new TH1D(Form("nume_phi%d",i),";#phi;Events",Nphiarr-1,phiarr);
    deno_phi[i] = new TH1D(Form("deno_phi%d",i),";#phi;Events",Nphiarr-1,phiarr);*/
  };
  
  //Int_t nevent = fcha->GetEntries();
  Int_t nevent = 100000;

  for(int i=0; i<nevent; i++){
    fcha->GetEntry(i);
    if(i%20000==0){cout<<">>>>> EVENT "<<i<<" / "<<fcha->GetEntries()<<" ("<<(int)(100.*i/fcha->GetEntries())<<"%)"<<endl;}
    
    for(int j=0; j<Reco_mu_size; j++){
      //Just getting the 4 momentum
      TLorentzVector *recomu4mom = (TLorentzVector*) Reco_mu_4mom->At(j);
      //SoftMuon Cuts
      Bool_t SingleCond = true;
      SingleCond = SingleCond && (Reco_mu_isGoodMuon[j]==1);
      SingleCond = SingleCond && (Reco_mu_nTrkWMea[j] > 5);
      SingleCond = SingleCond && (Reco_mu_nPixWMea[j] > 0);
      SingleCond = SingleCond && (fabs(Reco_mu_dxy[j]) < 0.3);
      SingleCond = SingleCond && (fabs(Reco_mu_dz[j]) < 20.);

      for(int k=kStart; k<Ntrig; k++){
        if( SingleCond//soft muon cut                
            && AcceptanceCut(recomu4mom)//Acceptance cut
          )
        {
          deno_p[k]->Fill(recomu4mom->Pt());
          double ptCut = getPtCut(trigname[k].c_str());
            if(recomu4mom->Pt()>ptCut) deno_e[k]->Fill(recomu4mom->Eta());//denominator for SingleMu trigs

          //Fill numerator if it passes an additional requirement
          if( ((Reco_mu_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) && ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))) 
          {
            nume_p[k]->Fill(recomu4mom->Pt());
            if(recomu4mom->Pt()>ptCut) nume_e[k]->Fill(recomu4mom->Eta());
          }//end of numerator fill
        }//end of SingleCond and Acceptance Cut
      }//end of k loop
    }//for Fill
  }//for nevent

  TGraphAsymmErrors *grp[10];
  TGraphAsymmErrors *gre[10];
  for(int k=kStart;k<Ntrig;k++){
    nume_p[k]->Divide(nume_p[k],deno_p[k],1,1,"B");
    nume_e[k]->Divide(nume_e[k],deno_e[k],1,1,"B");
    grp[k] = new TGraphAsymmErrors(nume_p[k]);
    gre[k] = new TGraphAsymmErrors(nume_e[k]);

    if(k-kStart<4){
      SetGraphStyle(grp[k],k-kStart,k-kStart);
      SetGraphStyle(gre[k],k-kStart,k-kStart);
    }
    else{
      SetGraphStyleOpen(grp[k],k-kStart-4,k-kStart-4,k-kStart-4);
      SetGraphStyleOpen(gre[k],k-kStart-4,k-kStart-4,k-kStart-4);
    }
  }
  for (int k=kStart;k<Ntrig;k++){
    grp[k]->GetYaxis()->SetRangeUser(0,1.3);
    grp[k]->GetXaxis()->SetLimits(ptmin,ptmax);
    grp[k]->GetXaxis()->SetRangeUser(ptmin,ptmax);
    grp[k]->GetXaxis()->SetTitle("p^{#mu}_{T} (GeV/c))");
    grp[k]->GetYaxis()->SetTitle("Efficiency");
    grp[k]->GetYaxis()->CenterTitle();
    grp[k]->GetXaxis()->CenterTitle();
    gre[k]->GetYaxis()->SetRangeUser(0,1.3);
    gre[k]->GetXaxis()->SetRangeUser(etamin,etamax);
    gre[k]->GetXaxis()->SetLimits(etamin,etamax);
    gre[k]->GetXaxis()->SetTitle("#eta^{#mu}");
    gre[k]->GetYaxis()->SetTitle("Efficiency");
    gre[k]->GetYaxis()->CenterTitle();
    gre[k]->GetXaxis()->CenterTitle();
  }

  //Save the Graphs in a root file
  TFile* outFile = new TFile("HLT_SingleMu_muGun_all.root","recreate");
  for(int k=kStart;k<Ntrig;k++){
    TString histnamep = Form("%s%s",trigname[k].c_str(),"_pt");
    TString histnamee = Form("%s%s",trigname[k].c_str(),"_eta");
    grp[k]->Write(histnamep);
    gre[k]->Write(histnamee);
  }
  outFile->Close();
}
