#define GluonHistosFill_cxx

#include "GluonHistosFill.h"

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

#include <vector>

void GluonHistosFill::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GluonHistosFill.C
//      root> GluonHistosFill t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:

fChain->SetBranchStatus("*",0);  // disable all branches
fChain->SetBranchStatus("*jetPt",1);  // activate branchname
fChain->SetBranchStatus("*genJetPt",1);
fChain->SetBranchStatus("*nGenJetPF",1);
fChain->SetBranchStatus("*isPhysG",1);
fChain->SetBranchStatus("*isPhysUDS",1);
fChain->SetBranchStatus("*jetEta",1);
fChain->SetBranchStatus("*jetPtOrder",1);
fChain->SetBranchStatus("*jetGirth",1);
fChain->SetBranchStatus("*jetMass",1);
fChain->SetBranchStatus("*jetGirth",1);
fChain->SetBranchStatus("*PF_dR",1);
fChain->SetBranchStatus("*PF_fromPV",1);
fChain->SetBranchStatus("*jetArea",1);

// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   const int NBINS = 62;
   double ptrange[NBINS+1] = {
      1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1588, 1784, 2000, 2238, 2500, 3500, 3832, 4252, 4713, 5220, 5777, 6389, 7000
   };

   const  int NPFBINS = 250;
   double pfrange[NPFBINS+1];
   for (int i = 0; i < NPFBINS+1; ++i) {
      pfrange[ i ] = i;
   }

   const int NMASSBINS = 200;
   double massrange[NMASSBINS+1];
   for (int i = 0; i < NMASSBINS+1; ++i) {
      massrange[ i ] = i;
   }

   #ifdef RECREATE_WEIGHTS
   FillWeightHistos(NBINS, ptrange, NPFBINS, pfrange, NMASSBINS, massrange);
   #endif

   #ifdef SINGLE_TREE
   TFile* f = TFile::Open("weightHistos_single.root");
   #else
   TFile* f = TFile::Open("weightHistos.root");
   #endif
   TH2D* gHist2;
   TH2D* qHist2;
   TH1D* nUEPF_perA_hist;
   TH2D* gHist2_sig;
   TH2D* qHist2_sig;

   f->GetObject("gluon_nGenPF_probs", gHist2);
   f->GetObject("quark_nGenPF_probs", qHist2);
   f->GetObject("UE_nPF_perA", nUEPF_perA_hist);
   f->GetObject("gluon_nGenPFSig_probs", gHist2_sig);
   f->GetObject("quark_nGenPFSig_probs", qHist2_sig);

   #ifdef SINGLE_TREE
   TFile* file = new TFile("outputGluonHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("outputGluonHistos.root", "recreate");
   #endif

   //All
   TProfile* PtResp = new TProfile("pt_resp", "", NBINS, ptrange);
   TProfile* NGenPF = new TProfile("nGenPF_prof", "", NBINS, ptrange);
   TProfile* NGenPF_wg = new TProfile("nGenPF_prof_wg", "", NBINS, ptrange);
   TProfile* PtRespNGenPF_wg = new TProfile("pt_resp_nGenJetPF_wg", "", NBINS, ptrange);
   TProfile* PtRespNGenPFSig_wg = new TProfile("pt_resp_nGenJetPFSig_wg", "", NBINS, ptrange);
   TProfile* NGenPF_wq = new TProfile("nGenPF_prof_wq", "", NBINS, ptrange);
   TProfile* PtRespNGenPF_wq = new TProfile("pt_resp_nGenJetPF_wq", "", NBINS, ptrange);
   TProfile* PtRespNGenPFSig_wq = new TProfile("pt_resp_nGenJetPFSig_wq", "", NBINS, ptrange);

   TH1D* PerPtBin = new TH1D("per_pt_bin", "", NBINS, ptrange);
   TH2D* PtResp_hist = new TH2D("pt_resp_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* PtRespNGenPF_wg_hist = new TH2D("pt_resp_nGenJetPF_wg_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* PtRespNGenPFSig_wg_hist = new TH2D("pt_resp_nGenJetPFSig_wg_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* PtRespNGenPF_wq_hist = new TH2D("pt_resp_nGenJetPF_wq_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* PtRespNGenPFSig_wq_hist = new TH2D("pt_resp_nGenJetPFSig_wq_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);


   // Gluons
   TProfile* gPtResp = new TProfile("gluon_pt_resp", "", NBINS, ptrange);
   TProfile* gNGenPF = new TProfile("gluon_nGenPF_prof", "", NBINS, ptrange);
   TProfile* gNGenPF_w = new TProfile("gluon_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* gPtRespNGenPF_w = new TProfile("gluon_pt_resp_nGenJetPF_w", "", NBINS, ptrange);
   TProfile* gPtRespNGenPFSig_w = new TProfile("gluon_pt_resp_nGenJetPFSig_w", "", NBINS, ptrange);

   TH2D* gNGenPF_w_hist = new TH2D("gluon_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* gPerPtBin = new TH1D("gluons_per_pt_bin", "", NBINS, ptrange);
   TH2D* gPtResp_hist = new TH2D("gluon_pt_resp_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* gPtResp_hist_norm = new TH2D("gluon_pt_resp_hist_norm", "", NBINS, ptrange, 1000, 0.5,1.5);
   TH2D* gPtRespNGenPF_w_hist = new TH2D("gluon_pt_resp_nGenJetPF_w_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* gPtRespNGenPFSig_w_hist = new TH2D("gluon_pt_resp_nGenJetPFSig_w_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);

   // Quarks
   TProfile* qPtResp = new TProfile("quark_pt_resp", "", NBINS, ptrange);
   TProfile* qNGenPF = new TProfile("quark_nGenPF_prof", "", NBINS, ptrange);
   TProfile* qNGenPF_w = new TProfile("quark_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* qPtRespNGenPF_w = new TProfile("quark_pt_resp_nGenJetPF_w", "", NBINS, ptrange);
   TProfile* qPtRespNGenPFSig_w = new TProfile("quark_pt_resp_nGenJetPFSig_w", "", NBINS, ptrange);

   TH2D* qNGenPF_w_hist = new TH2D("quark_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* qPerPtBin = new TH1D("quarks_per_pt_bin", "", NBINS, ptrange);
   TH2D* qPtResp_hist = new TH2D("quark_pt_resp_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* qPtResp_hist_norm = new TH2D("quark_pt_resp_hist_norm", "", NBINS, ptrange, 1000, 0.5,1.5);
   TH2D* qPtRespNGenPF_w_hist = new TH2D("quark_pt_resp_nGenJetPF_w_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);
   TH2D* qPtRespNGenPFSig_w_hist = new TH2D("quark_pt_resp_nGenJetPFSig_w_hist", "", NBINS, ptrange, 1000, 0.5, 1.5);

   TH2D* nGenJetPF_probsRatio = new TH2D("nGenJetPF_likelyhood", "", NBINS, ptrange, NPFBINS, pfrange);

   TH2D* nGenPF_probs_sum = (TH2D*)gHist2->Clone("nGenPF_probs_sum");

   nGenPF_probs_sum->Add(qHist2);
   nGenJetPF_probsRatio->Divide(qHist2, nGenPF_probs_sum);

   for (int xb = 1; xb != nGenJetPF_probsRatio->GetNbinsX()+1; ++xb) {
      for(int yb = 1; yb != nGenJetPF_probsRatio->GetNbinsY()+1; ++yb) {
            nGenJetPF_probsRatio->SetBinContent(xb, yb, (nGenJetPF_probsRatio->GetBinContent(xb, yb)-0.5)*2);
      }
   }

   double resp;
   double w;
   double gprob;
   double qprob;

   double nUEPF_per_A = nUEPF_perA_hist->GetMean();
   int nPFsum;
   double dR_in = 0.8;
   double dR_out = 1.0;
   double nUEPF;
   double nGenPFSig;

   double w_sig;
   double gprob_sig;
   double qprob_sig;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
         if (recoPlotted) {
            usedPlotPt = jetPt;
         }
         else {
            usedPlotPt = genJetPt;
         }
         if (recoWeighted) {
            usedWeightPt = jetPt;
         }
         else {
            usedWeightPt = genJetPt;
         }

         resp = jetPt/genJetPt;
         gprob = gHist2->GetBinContent(gHist2->FindBin(usedWeightPt, nGenJetPF));
         qprob = qHist2->GetBinContent(qHist2->FindBin(usedWeightPt, nGenJetPF));

         if (isPFSig_perEvent) {
               nPFsum = 0;
               for (int i = 0; i != nPF; ++i) {
               if (PF_dR[i] > dR_in && PF_dR[i] < dR_out && PF_fromPV[i] == 3) {
                  ++nPFsum;
               }
            }
            nUEPF_per_A = 2.75 * nPFsum / (M_PI * (dR_out*dR_out - dR_in*dR_in));
         }
         nUEPF = jetArea*nUEPF_per_A;
         nGenPFSig = nGenJetPF - nUEPF;
         gprob_sig = gHist2_sig->GetBinContent(gHist2_sig->FindBin(usedWeightPt, nGenPFSig));
         qprob_sig = qHist2_sig->GetBinContent(qHist2_sig->FindBin(usedWeightPt, nGenPFSig));
         if (isPhysG) {
            gPtResp->Fill(usedPlotPt, resp);
            gPtResp_hist->Fill(usedPlotPt, resp);
            gNGenPF->Fill(usedPlotPt, nGenJetPF);
            gPerPtBin->Fill(usedPlotPt);
            if (gprob_sig > 0 && qprob_sig > 0) {
               w_sig = qprob_sig/gprob_sig;
               gPtRespNGenPFSig_w->Fill(usedPlotPt, resp, w_sig);
               gPtRespNGenPFSig_w_hist->Fill(usedPlotPt, resp, w_sig);
            }
            if (gprob > 0 && qprob > 0) {
               w = qprob/gprob;
               gPtRespNGenPF_w->Fill(usedPlotPt, resp, w);
               gPtRespNGenPF_w_hist->Fill(usedPlotPt, resp, w);
            }
            else w = 0;
            gNGenPF_w_hist->Fill(usedPlotPt, nGenJetPF, w);
            gNGenPF_w->Fill(usedPlotPt, nGenJetPF, w);
         }
         else {
            if (isPhysUDS) {
               qPtResp->Fill(usedPlotPt, resp);
               qPtResp_hist->Fill(usedPlotPt, resp);
               qNGenPF->Fill(usedPlotPt, nGenJetPF);
               qPerPtBin->Fill(usedPlotPt);
               if (gprob_sig > 0 && qprob_sig > 0) {
                  w_sig = gprob_sig/qprob_sig;
                  qPtRespNGenPFSig_w->Fill(usedPlotPt, resp, w_sig);
                  qPtRespNGenPFSig_w_hist->Fill(usedPlotPt, resp, w_sig);
               }
               if (gprob > 0 && qprob > 0) {
                  w = gprob/qprob;
                  qPtRespNGenPF_w->Fill(usedPlotPt, resp, w);
                  qPtRespNGenPF_w_hist->Fill(usedPlotPt, resp, w);
               }
               else w = 0;
               qNGenPF_w_hist->Fill(usedPlotPt, nGenJetPF, w);
               qNGenPF_w->Fill(usedPlotPt, nGenJetPF, w);
            }
         }
         //All
         PtResp->Fill(usedPlotPt, resp);
         PtResp_hist->Fill(usedPlotPt, resp);
         if (gprob > 0 && qprob > 0) {
            PtRespNGenPF_wg->Fill(usedPlotPt, resp, qprob/gprob);
            PtRespNGenPF_wg_hist->Fill(usedPlotPt, resp, qprob/gprob);

            PtRespNGenPF_wq->Fill(usedPlotPt, resp, gprob/qprob);
            PtRespNGenPF_wq_hist->Fill(usedPlotPt, resp, gprob/qprob);
         }
         if (gprob_sig > 0 && qprob_sig > 0) {
            PtRespNGenPFSig_wg->Fill(usedPlotPt, resp, qprob_sig/gprob_sig);
            PtRespNGenPFSig_wg_hist->Fill(usedPlotPt, resp, qprob_sig/gprob_sig);

            PtRespNGenPFSig_wq->Fill(usedPlotPt, resp, gprob_sig/qprob_sig);
            PtRespNGenPFSig_wq_hist->Fill(usedPlotPt, resp, gprob_sig/qprob_sig);
         }

         PerPtBin->Fill(usedPlotPt);
   }
}

double fit_up;
double fit_down;

PtResp_hist->SetAxisRange(fit_down, fit_up, "Y");
PtRespNGenPF_wg_hist->SetAxisRange(fit_down, fit_up, "Y");
PtRespNGenPFSig_wg_hist->SetAxisRange(fit_down, fit_up, "Y");
PtRespNGenPF_wq_hist->SetAxisRange(fit_down, fit_up, "Y");
PtRespNGenPFSig_wq_hist->SetAxisRange(fit_down, fit_up, "Y");

gPtResp_hist->SetAxisRange(fit_down, fit_up, "Y");
gPtRespNGenPF_w_hist->SetAxisRange(fit_down, fit_up, "Y");
gPtRespNGenPFSig_w_hist->SetAxisRange(fit_down, fit_up, "Y");

qPtResp_hist->SetAxisRange(fit_down, fit_up, "Y");
qPtRespNGenPF_w_hist->SetAxisRange(fit_down, fit_up, "Y");
qPtRespNGenPFSig_w_hist->SetAxisRange(fit_down, fit_up, "Y");

PtResp_hist->FitSlicesY();
PtRespNGenPF_wg_hist->FitSlicesY();
PtRespNGenPFSig_wg_hist->FitSlicesY();
PtRespNGenPF_wq_hist->FitSlicesY();
PtRespNGenPFSig_wq_hist->FitSlicesY();

gPtResp_hist->FitSlicesY();
gPtRespNGenPF_w_hist->FitSlicesY();
gPtRespNGenPFSig_w_hist->FitSlicesY();

qPtResp_hist->FitSlicesY();
qPtRespNGenPF_w_hist->FitSlicesY();
qPtRespNGenPFSig_w_hist->FitSlicesY();

CalculateProbs(gPtResp_hist, qPtResp_hist, gPtResp_hist_norm, qPtResp_hist_norm);

file->Write();
file->Close();

f->Close();
}


void GluonHistosFill::FillWeightHistos(int nptbins, double* ptrange, int npfbins, double* pfrange, int nMassBins, double* massrange)
{
   #ifdef SINGLE_TREE
   TFile* file = new TFile("./weightHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("./weightHistos.root", "recreate");
   #endif

   //Gluon
   TH2D* gNGenPF = new TH2D("gluon_nGenPF_hist", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* gGenJetMass = new TH2D("gluon_genJetMass", "", nptbins, ptrange, nMassBins, massrange);
   TH2D* gJetGirth = new TH2D("gluon_jetGirth", "", nptbins, ptrange, 100, 0., 0.5);

   TH2D* gNGenPF_probs = new TH2D("gluon_nGenPF_probs", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* gGenJetMass_probs = new TH2D("gluon_genJetMass_probs", "", nptbins, ptrange, nMassBins, massrange);
   TH2D* gJetGirth_probs = new TH2D("gluons_jetGirth_probs", "", nptbins, ptrange, 100, 0., 0.5);
   TH2D* gNGenPFSig = new TH2D("gluon_nGenPFSig", "", nptbins, ptrange, 300, -50, 250);
   TH2D* gNGenPFSig_probs = new TH2D("gluon_nGenPFSig_probs", "", nptbins, ptrange, 300, -50, 250);
   TProfile* gNGenPFSig_prof = new TProfile("gluon_nGenPFSig_prof", "", nptbins, ptrange);

   //Quarks
   TH2D* qNGenPF = new TH2D("quark_nGenPF_hist", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* qGenJetMass = new TH2D("quark_genJetMass", "", nptbins, ptrange, nMassBins, pfrange);
   TH2D* qJetGirth = new TH2D("quark_jetGirth", "", nptbins, ptrange, 100, 0., 0.5);

   TH2D* qNGenPF_probs = new TH2D("quark_nGenPF_probs", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* qGenJetMass_probs = new TH2D("quark_genJetMass_probs", "", nptbins, ptrange, nMassBins, massrange);
   TH2D* qJetGirth_probs = new TH2D("quarks_jetGirth_probs", "", nptbins, ptrange, 100, 0., 0.5);
   TH2D* qNGenPFSig = new TH2D("quark_nGenPFSig", "", nptbins, ptrange, 300, -50, 250);
   TH2D* qNGenPFSig_probs = new TH2D("quark_nGenPFSig_probs", "", nptbins, ptrange, 300, -50, 250);
   TProfile* qNGenPFSig_prof = new TProfile("quark_nGenPFSig_prof", "", nptbins, ptrange);

   int nPFsum;
   double nUEPF_per_A;
   double nUEPF;
   double dR_in = 0.8;
   double dR_out = 1.0;
   double nGenPFSig;

   TH1D* nUEPF_perA_hist = new TH1D("UE_nPF_perA", "", int(80./2.75), 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_nGenPF_hist = new TH2D("nUEPF_perA_nGenPF_hist", "", npfbins, pfrange, int(80./2.75), 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_nGenPF_probs = new TH2D("nUEPF_perA_nGenPF_probs", "", npfbins, pfrange, int(80./2.75), 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_genJetPt_hist = new TH2D("nUEPF_perA_genJetPt_hist", "", nptbins, ptrange, int(80./2.75), 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_genJetPt_probs = new TH2D("nUEPF_perA_genJetPt_probs", "", nptbins, ptrange, int(80./2.75), 0, int(80./2.75)*2.75);
   TProfile* nUEPF_nGenPF = new TProfile("nUEPF_nGenPF", "", npfbins, pfrange);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
         if (recoWeighted) {
            usedWeightPt = jetPt;
         }
         else usedWeightPt = genJetPt;

         if (isPhysG) {
            gNGenPF->Fill(usedWeightPt, nGenJetPF);
            gGenJetMass->Fill(usedWeightPt, genJetMass);
            gJetGirth->Fill(usedWeightPt, jetGirth);
            nPFsum = 0;
            for (int i = 0; i != nPF; ++i) {
                  if (PF_dR[i] > dR_in && PF_dR[i] < dR_out && PF_fromPV[i] == 3) {
                     ++nPFsum;
                  }
            }
            nUEPF_per_A = 2.75 * nPFsum / (M_PI * (dR_out*dR_out - dR_in*dR_in));
            nUEPF_perA_hist->Fill(nUEPF_per_A);
            nUEPF_perA_genJetPt_hist->Fill(usedWeightPt, nUEPF_per_A);
            nUEPF_perA_nGenPF_hist->Fill(nGenJetPF, nUEPF_per_A);
            if (isPFSig_perEvent) {
               nUEPF = jetArea*nUEPF_per_A;
               nGenPFSig = nGenJetPF - nUEPF;
               gNGenPFSig->Fill(usedWeightPt, nGenPFSig);
            }
         }
         else {
            if (isPhysUDS) {
               qNGenPF->Fill(usedWeightPt, nGenJetPF);
               qGenJetMass->Fill(usedWeightPt, genJetMass);
               qJetGirth->Fill(usedWeightPt, jetGirth);
               nPFsum = 0;
               for (int i = 0; i != nPF; ++i) {
                  if (PF_dR[i] > dR_in && PF_dR[i] < dR_out && PF_fromPV[i] == 3) {
                     ++nPFsum;
                  }
               }
               nUEPF_per_A = 2.75 * nPFsum / (M_PI * (dR_out*dR_out - dR_in*dR_in));
               nUEPF_perA_hist->Fill(nUEPF_per_A);
               nUEPF_perA_genJetPt_hist->Fill(usedWeightPt, nUEPF_per_A);
               nUEPF_perA_nGenPF_hist->Fill(nGenJetPF, nUEPF_per_A);
               if (isPFSig_perEvent) {
                  nUEPF = jetArea*nUEPF_per_A;
                  nGenPFSig = nGenJetPF - nUEPF;
                  qNGenPFSig->Fill(usedWeightPt, nGenPFSig);
               }
            }
         }
      }
   }

   cout << "first loop done" << endl;

   nUEPF_per_A = nUEPF_perA_hist->GetMean();
   cout << nUEPF_per_A << endl;

   if (!isPFSig_perEvent) {
      nentries = fChain->GetEntriesFast();
      nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
                     if (recoWeighted) {
               usedWeightPt = jetPt;
            }
            else usedWeightPt = genJetPt;

            if (isPhysG) {
               nUEPF = jetArea*nUEPF_per_A;
               nGenPFSig = nGenJetPF - nUEPF;
               gNGenPFSig->Fill(usedWeightPt, nGenPFSig);
               gNGenPFSig_prof->Fill(usedWeightPt, nGenPFSig);
            }
            else {
               if (isPhysUDS) {
                  nUEPF = jetArea*nUEPF_per_A;
                  nGenPFSig = nGenJetPF - nUEPF;
                  qNGenPFSig->Fill(usedWeightPt, nGenPFSig);
                  qNGenPFSig_prof->Fill(usedWeightPt, nGenPFSig);
               }
            }
         }
      }
   }

   gNGenPF->Smooth(1, "k5b");
   gNGenPFSig->Smooth(1, "k5b");
   qNGenPF->Smooth(1, "k5b");
   qNGenPFSig->Smooth(1, "k5b");

   CalculateProbs(nUEPF_perA_genJetPt_hist, nullptr, nUEPF_perA_genJetPt_probs, nullptr);
   CalculateProbs(nUEPF_perA_nGenPF_hist, nullptr, nUEPF_perA_nGenPF_probs, nullptr);

   CalculateProbs(gNGenPF, qNGenPF, gNGenPF_probs, qNGenPF_probs);
   CalculateProbs(gNGenPFSig, qNGenPFSig, gNGenPFSig_probs, qNGenPFSig_probs);

   CalculateProbs(gGenJetMass, qGenJetMass, gGenJetMass_probs, qGenJetMass_probs);
   CalculateProbs(gJetGirth, qJetGirth, gJetGirth_probs, qJetGirth_probs);

   cout << "probability calculated" << endl;

   file->Write();
   file->Close();
}

void GluonHistosFill::CalculateProbs(TH2D* gHist, TH2D* qHist, TH2D* gProbs, TH2D* qProbs) 
{
   double gIntegral;
   double qIntegral;
   for (int xb = 1; xb != gHist->GetNbinsX()+1; ++xb) {
      gIntegral = gHist->Integral(xb, xb, 0, gHist->GetNbinsY());
      if (qHist != nullptr) {
         qIntegral = qHist->Integral(xb, xb, 0, qHist->GetNbinsY());
      }
      if (gIntegral > 0 || (qIntegral > 0 && qProbs != nullptr)) {
         for (int yb = 1; yb != gHist->GetNbinsY()+1; ++yb) {
            if (gIntegral > 0) {
               gProbs->SetBinContent(xb, yb, gHist->GetBinContent(xb, yb) / gIntegral);
            }
            if (qIntegral > 0 && qHist != nullptr) {
               qProbs->SetBinContent(xb, yb, qHist->GetBinContent(xb, yb) / qIntegral);
            }
         }
      }
   }
}