//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 28 12:03:51 2021 by ROOT version 6.24/00
// from TChain AK4jets/jetTree/
//////////////////////////////////////////////////////////

#ifndef GluonHistosFill_h
#define GluonHistosFill_h

#define SINGLE_TREE
//#define RECREATE_WEIGHTS

#include <iostream>

#include <TH2.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

// Header file for the classes stored in the TTree if any.

const int N_PF = 1000;
const int N_GenJetPF = 500;


class GluonHistosFill {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         jetPt;
   Float_t         jetEta;
   Float_t         jetPhi;
   Float_t         jetMass;
   Float_t         jetGirth;
   Float_t         jetArea;
   Float_t         jetRawPt;
   Float_t         jetRawMass;
   Int_t           jetLooseID;
   Int_t           jetTightID;
   Int_t           jetGenMatch;
   Float_t         jetQGl;
   Float_t         QG_ptD;
   Float_t         QG_axis2;
   Int_t           QG_mult;
   Int_t           partonFlav;
   Int_t           hadronFlav;
   Int_t           physFlav;
   Int_t           isPhysUDS;
   Int_t           isPhysG;
   Int_t           isPhysOther;
   Int_t           isPartonUDS;
   Int_t           isPartonG;
   Int_t           isPartonOther;
   Int_t           jetChargedHadronMult;
   Int_t           jetNeutralHadronMult;
   Int_t           jetChargedMult;
   Int_t           jetNeutralMult;
   Int_t           jetMult;
   Int_t           nPF;
   Float_t         PF_pT[N_PF];   //[nPF]
   Float_t         PF_dR[N_PF];   //[nPF]
   Float_t         PF_dTheta[N_PF];   //[nPF]
   Float_t         PF_dPhi[N_PF];   //[nPF]
   Float_t         PF_dEta[N_PF];   //[nPF]
   Float_t         PF_mass[N_PF];   //[nPF]
   Int_t           PF_id[N_PF];   //[nPF]
   Int_t           PF_fromPV[N_PF];   //[nPF]
   Int_t           PF_fromAK4Jet[N_PF];   //[nPF]
   Float_t         genJetPt;
   Float_t         genJetEta;
   Float_t         genJetPhi;
   Float_t         genJetMass;
   Int_t           nGenJetPF;
   Float_t         genJetPF_pT[N_GenJetPF];   //[nGenJetPF]
   Float_t         genJetPF_dR[N_GenJetPF];   //[nGenJetPF]
   Float_t         genJetPF_dTheta[N_GenJetPF];   //[nGenJetPF]
   Float_t         genJetPF_mass[N_GenJetPF];   //[nGenJetPF]
   Int_t           genJetPF_id[N_GenJetPF];   //[nGenJetPF]
   Int_t           eventJetMult;
   Int_t           jetPtOrder;
   Float_t         dPhiJetsLO;
   Float_t         dEtaJetsLO;
   Float_t         alpha;
   ULong64_t       event;
   Int_t           run;
   Int_t           lumi;
   Float_t         pthat;
   Float_t         eventWeight;
   Float_t         rhoAll;
   Float_t         rhoCentral;
   Float_t         rhoCentralNeutral;
   Float_t         rhoCentralChargedPileUp;
   Int_t           PV_npvsGood;
   Int_t           Pileup_nPU;
   Float_t         Pileup_nTrueInt;

   // List of branches
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetGirth;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawMass;   //!
   TBranch        *b_jetLooseID;   //!
   TBranch        *b_jetTightID;   //!
   TBranch        *b_jetGenMatch;   //!
   TBranch        *b_jetQGl;   //!
   TBranch        *b_QG_ptD;   //!
   TBranch        *b_QG_axis2;   //!
   TBranch        *b_QG_mult;   //!
   TBranch        *b_partonFlav;   //!
   TBranch        *b_hadronFlav;   //!
   TBranch        *b_physFlav;   //!
   TBranch        *b_isPhysUDS;   //!
   TBranch        *b_isPhysG;   //!
   TBranch        *b_isPhysOther;   //!
   TBranch        *b_isPartonUDS;   //!
   TBranch        *b_isPartonG;   //!
   TBranch        *b_isPartonOther;   //!
   TBranch        *b_jetChargedHadronMult;   //!
   TBranch        *b_jetNeutralHadronMult;   //!
   TBranch        *b_jetChargedMult;   //!
   TBranch        *b_jetNeutralMult;   //!
   TBranch        *b_jetMult;   //!
   TBranch        *b_nPF;   //!
   TBranch        *b_PF_pT;   //!
   TBranch        *b_PF_dR;   //!
   TBranch        *b_PF_dTheta;   //!
   TBranch        *b_PF_dPhi;   //!
   TBranch        *b_PF_dEta;   //!
   TBranch        *b_PF_mass;   //!
   TBranch        *b_PF_id;   //!
   TBranch        *b_PF_fromPV;   //!
   TBranch        *b_PF_fromAK4Jet;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetMass;   //!
   TBranch        *b_nGenJetPF;   //!
   TBranch        *b_genJetPF_pT;   //!
   TBranch        *b_genJetPF_dR;   //!
   TBranch        *b_genJetPF_dTheta;   //!
   TBranch        *b_genJetPF_mass;   //!
   TBranch        *b_genJetPF_id;   //!
   TBranch        *b_eventJetMult;   //!
   TBranch        *b_jetPtOrder;   //!
   TBranch        *b_dPhiJetsLO;   //!
   TBranch        *b_dEtaJetsLO;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_eventWeight;   //!
   TBranch        *b_rhoAll;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_rhoCentralNeutral;   //!
   TBranch        *b_rhoCentralChargedPileUp;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_nTrueInt;   //!

   GluonHistosFill(TTree *tree=0);
   virtual ~GluonHistosFill();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void             FillWeightHistos(int nptbins, double* ptrange, int npfbins, double* pfrange, int nMassBins, double* massrange);
   void             CalculateProbs(TH2D* gHist, TH2D* qHist, TH2D* gProbs, TH2D* qProbs);
};

#endif

#ifdef GluonHistosFill_cxx
GluonHistosFill::GluonHistosFill(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_415.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_415.root");
      }
      f->GetObject("AK4jets/jetTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("AK4jets/jetTree","");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_1.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_2.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_3.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_4.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_5.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_6.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_7.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_8.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_9.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_10.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_11.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_12.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_13.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_14.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_15.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_16.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_17.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_18.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_19.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_20.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_21.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_22.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_23.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_24.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_25.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_26.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_27.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_28.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_29.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_30.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_31.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_32.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_33.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_34.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_35.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_36.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_37.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_38.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_39.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_40.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_41.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_42.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_43.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_44.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_45.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_46.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_47.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_48.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_49.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_50.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_51.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_52.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_53.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_54.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_55.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_56.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_57.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_58.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_59.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_60.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_61.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_62.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_63.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_64.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_65.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_66.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_67.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_68.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_69.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_70.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_71.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_72.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_73.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_74.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_75.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_76.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_77.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_78.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_79.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_80.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_81.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_82.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_83.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_84.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_85.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_86.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_87.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_88.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_89.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_90.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_91.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_92.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_93.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_94.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_95.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_96.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_97.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_98.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_99.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_100.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_101.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_102.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_103.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_104.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_105.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_106.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_107.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_108.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_109.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_110.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_111.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_112.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_113.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_114.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_115.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_116.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_117.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_118.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_119.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_120.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_121.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_122.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_123.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_124.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_125.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_126.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_127.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_128.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_129.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_130.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_131.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_132.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_133.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_134.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_135.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_136.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_137.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_138.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_139.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_140.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_141.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_142.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_143.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_144.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_145.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_146.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_147.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_148.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_149.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_150.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_151.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_152.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_153.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_154.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_155.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_156.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_157.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_158.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_159.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_160.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_161.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_162.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_163.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_164.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_165.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_166.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_167.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_168.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_169.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_170.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_171.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_172.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_173.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_174.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_175.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_176.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_177.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_178.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_179.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_180.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_181.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_182.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_183.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_184.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_185.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_186.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_187.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_188.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_189.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_190.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_191.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_192.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_193.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_194.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_195.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_196.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_197.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_198.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_199.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_200.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_201.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_202.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_203.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_204.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_205.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_206.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_207.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_208.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_209.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_210.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_211.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_212.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_213.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_214.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_215.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_216.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_217.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_218.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_219.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_220.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_221.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_222.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_223.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_224.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_225.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_226.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_227.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_228.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_229.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_230.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_231.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_232.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_233.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_234.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_235.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_236.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_237.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_238.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_239.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_240.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_241.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_242.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_243.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_244.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_245.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_246.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_247.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_248.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_249.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_250.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_251.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_252.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_253.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_254.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_255.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_256.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_257.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_258.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_259.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_260.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_261.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_262.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_263.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_264.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_265.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_266.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_267.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_268.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_269.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_270.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_271.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_272.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_273.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_274.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_275.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_276.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_277.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_278.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_279.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_280.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_281.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_282.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_283.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_284.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_285.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_286.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_287.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_288.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_289.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_290.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_291.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_292.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_293.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_294.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_295.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_296.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_297.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_298.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_299.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_300.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_301.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_302.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_303.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_304.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_305.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_306.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_307.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_308.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_309.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_310.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_311.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_312.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_313.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_314.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_315.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_316.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_317.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_318.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_319.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_320.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_321.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_322.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_323.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_324.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_325.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_326.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_327.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_328.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_329.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_330.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_331.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_332.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_333.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_334.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_335.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_336.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_337.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_338.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_339.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_340.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_341.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_342.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_343.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_344.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_345.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_346.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_347.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_348.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_349.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_350.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_351.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_352.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_353.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_354.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_355.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_356.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_357.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_358.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_359.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_360.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_361.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_362.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_363.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_364.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_365.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_366.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_367.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_368.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_369.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_370.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_371.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_372.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_373.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_374.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_375.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_376.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_377.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_378.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_379.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_380.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_381.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_382.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_383.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_384.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_385.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_386.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_387.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_388.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_389.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_390.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_391.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_392.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_393.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_394.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_395.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_396.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_397.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_398.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_399.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_400.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_401.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_402.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_403.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_404.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_405.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_406.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_407.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_408.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_409.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_410.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_411.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_412.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_413.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_414.root/AK4jets/jetTree");
      chain->Add("/work/gluon_tuples/multicrab_ConfFile_cfg_v10625_2018MC_20210526T0830/QCD_Pt_15to7000/results/JetNtuple_415.root/AK4jets/jetTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

GluonHistosFill::~GluonHistosFill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GluonHistosFill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GluonHistosFill::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GluonHistosFill::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetGirth", &jetGirth, &b_jetGirth);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawMass", &jetRawMass, &b_jetRawMass);
   fChain->SetBranchAddress("jetLooseID", &jetLooseID, &b_jetLooseID);
   fChain->SetBranchAddress("jetTightID", &jetTightID, &b_jetTightID);
   fChain->SetBranchAddress("jetGenMatch", &jetGenMatch, &b_jetGenMatch);
   fChain->SetBranchAddress("jetQGl", &jetQGl, &b_jetQGl);
   fChain->SetBranchAddress("QG_ptD", &QG_ptD, &b_QG_ptD);
   fChain->SetBranchAddress("QG_axis2", &QG_axis2, &b_QG_axis2);
   fChain->SetBranchAddress("QG_mult", &QG_mult, &b_QG_mult);
   fChain->SetBranchAddress("partonFlav", &partonFlav, &b_partonFlav);
   fChain->SetBranchAddress("hadronFlav", &hadronFlav, &b_hadronFlav);
   fChain->SetBranchAddress("physFlav", &physFlav, &b_physFlav);
   fChain->SetBranchAddress("isPhysUDS", &isPhysUDS, &b_isPhysUDS);
   fChain->SetBranchAddress("isPhysG", &isPhysG, &b_isPhysG);
   fChain->SetBranchAddress("isPhysOther", &isPhysOther, &b_isPhysOther);
   fChain->SetBranchAddress("isPartonUDS", &isPartonUDS, &b_isPartonUDS);
   fChain->SetBranchAddress("isPartonG", &isPartonG, &b_isPartonG);
   fChain->SetBranchAddress("isPartonOther", &isPartonOther, &b_isPartonOther);
   fChain->SetBranchAddress("jetChargedHadronMult", &jetChargedHadronMult, &b_jetChargedHadronMult);
   fChain->SetBranchAddress("jetNeutralHadronMult", &jetNeutralHadronMult, &b_jetNeutralHadronMult);
   fChain->SetBranchAddress("jetChargedMult", &jetChargedMult, &b_jetChargedMult);
   fChain->SetBranchAddress("jetNeutralMult", &jetNeutralMult, &b_jetNeutralMult);
   fChain->SetBranchAddress("jetMult", &jetMult, &b_jetMult);
   fChain->SetBranchAddress("nPF", &nPF, &b_nPF);
   fChain->SetBranchAddress("PF_pT", PF_pT, &b_PF_pT);
   fChain->SetBranchAddress("PF_dR", PF_dR, &b_PF_dR);
   fChain->SetBranchAddress("PF_dTheta", PF_dTheta, &b_PF_dTheta);
   fChain->SetBranchAddress("PF_dPhi", PF_dPhi, &b_PF_dPhi);
   fChain->SetBranchAddress("PF_dEta", PF_dEta, &b_PF_dEta);
   fChain->SetBranchAddress("PF_mass", PF_mass, &b_PF_mass);
   fChain->SetBranchAddress("PF_id", PF_id, &b_PF_id);
   fChain->SetBranchAddress("PF_fromPV", PF_fromPV, &b_PF_fromPV);
   fChain->SetBranchAddress("PF_fromAK4Jet", PF_fromAK4Jet, &b_PF_fromAK4Jet);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetMass", &genJetMass, &b_genJetMass);
   fChain->SetBranchAddress("nGenJetPF", &nGenJetPF, &b_nGenJetPF);
   fChain->SetBranchAddress("genJetPF_pT", genJetPF_pT, &b_genJetPF_pT);
   fChain->SetBranchAddress("genJetPF_dR", genJetPF_dR, &b_genJetPF_dR);
   fChain->SetBranchAddress("genJetPF_dTheta", genJetPF_dTheta, &b_genJetPF_dTheta);
   fChain->SetBranchAddress("genJetPF_mass", genJetPF_mass, &b_genJetPF_mass);
   fChain->SetBranchAddress("genJetPF_id", genJetPF_id, &b_genJetPF_id);
   fChain->SetBranchAddress("eventJetMult", &eventJetMult, &b_eventJetMult);
   fChain->SetBranchAddress("jetPtOrder", &jetPtOrder, &b_jetPtOrder);
   fChain->SetBranchAddress("dPhiJetsLO", &dPhiJetsLO, &b_dPhiJetsLO);
   fChain->SetBranchAddress("dEtaJetsLO", &dEtaJetsLO, &b_dEtaJetsLO);
   fChain->SetBranchAddress("alpha", &alpha, &b_alpha);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
   fChain->SetBranchAddress("rhoAll", &rhoAll, &b_rhoAll);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("rhoCentralNeutral", &rhoCentralNeutral, &b_rhoCentralNeutral);
   fChain->SetBranchAddress("rhoCentralChargedPileUp", &rhoCentralChargedPileUp, &b_rhoCentralChargedPileUp);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   Notify();
}

Bool_t GluonHistosFill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GluonHistosFill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GluonHistosFill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef GluonHistosFill_cxx