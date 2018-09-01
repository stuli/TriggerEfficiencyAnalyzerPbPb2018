#ifndef SetHLT_C
#define SetHLT_C

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>

struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };

struct TriggerList { int HLT_HIL1DoubleMu0_v1, HLT_HIL2DoubleMu0_v1, HLT_HIL3DoubleMu0_v2, HLT_L1DoubleMu0_v1, HLT_L1SingleMuOpen_v1, HLT_L1SingleMu3_v1, HLT_L1SingleMu5_v1, HLT_L1SingleMu7_v1, HLT_HIL2Mu7_v1, HLT_HIL3Mu3_v1, HLT_HIL3Mu5_v1, HLT_HIL3Mu7_v2; };
struct TriggerBranch { TBranch *b_HLT_HIL1DoubleMu0_v1, *b_HLT_HIL2DoubleMu0_v1, *b_HLT_HIL3DoubleMu0_v2, *b_HLT_L1DoubleMu0_v1, *b_HLT_L1SingleMuOpen_v1, *b_HLT_L1SingleMu3_v1, *b_HLT_L1SingleMu5_v1, *b_HLT_L1SingleMu7_v1, *b_HLT_HIL2Mu7_v1, *b_HLT_HIL3Mu3_v1, *b_HLT_HIL3Mu5_v1, *b_HLT_HIL3Mu7_v2; };

int kPromptJpsiLowPt = 0;
int kPromptJpsiHighPt = 1;
int kNonPromptJpsi = 2;
int kUpsione = 3;
int kSingleMu = 4;
int kZboson = 5;
int kXeXe=6;
  
const int Ntrig = 11;
/*
vector<string> trigname={
  "HLT_HIL1DoubleMuOpen_v1",
  "HLT_HIL1DoubleMuOpen_OS_v1",
  "HLT_HIL1DoubleMuOpen_SS_v1",
  "HLT_HIL1DoubleMuOpen_BptxAND_v1",
  "HLT_HIL1DoubleMuOpen_OS_BptxAND_v1",
  "HLT_HIL1DoubleMuOpen_SS_BptxAND_v1",
  "HLT_HIL1DoubleMuOpen_Cent10100_BptxAND_v1",
  "HLT_HIL1DoubleMuOpen_Cent50100_BptxAND_v1",
  "HLT_HIL1DoubleMu0_v1",
  "HLT_HIL1DoubleMu0_BptxAND_v1",
  "HLT_HIL1DoubleMu0_Cent10100_BptxAND_v1",
  "HLT_HIL1DoubleMu0_Cent30100_BptxAND_v1",
  "HLT_HIL1DoubleMu0_Cent50100_BptxAND_v1",
  "HLT_HIL2DoubleMu0_v1",
  "HLT_HIL3DoubleMu0_v2",
  "HLT_L1DoubleMu0_v1",
  "HLT_L1SingleMuOpen_v1",
  "HLT_L1SingleMu3_v1",
  "HLT_L1SingleMu5_v1",
  "HLT_L1SingleMu7_v1",
  "HLT_HIL2Mu7_v1",
  "HLT_HIL3Mu3_v1",
  "HLT_HIL3Mu5_v1",
  "HLT_HIL3Mu7_v2",
  "HLT_HIL3Mu12_v1",
  "HLT_HIL3Mu15_v1",
  "HLT_HIL3Mu20_v1"
};
*/
vector<string> trigname={
  "HLT_HIL1DoubleMu0_v1",
  "HLT_HIL2DoubleMu0_v1",
  "HLT_HIL3DoubleMu0_v2",
  "HLT_L1DoubleMu0_v1",
  "HLT_L1SingleMuOpen_v1",
  "HLT_L1SingleMu3_v1",
  "HLT_L1SingleMu5_v1",
  "HLT_L1SingleMu7_v1",
  "HLT_HIL2Mu7_v1",
  "HLT_HIL3Mu3_v1",
  "HLT_HIL3Mu5_v1",
  "HLT_HIL3Mu7_v2",
};

TString getSampleID( int sampleID ) {
  if ( sampleID == kPromptJpsiLowPt ) return "PR_Jpsi_lowPt";
  else if ( sampleID == kPromptJpsiHighPt ) return "PR_Jpsi_highPt";
  else if ( sampleID == kNonPromptJpsi ) return "NP_Jpsi";
  else if ( sampleID == kUpsione) return "Upsione";
  else if ( sampleID == kSingleMu) return "SingleMuGun";
  else if ( sampleID == kZboson) return "Z_Boson";
  else if ( sampleID == kXeXe) return "Z_Boson";
  else return "none";
}


#endif


