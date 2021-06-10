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
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

/*
   const int NBINS = 79;
   double ptrange[NBINS+1] = {
      1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
   };
*/

   const int NBINS = 62;
   double ptrange[NBINS+1] = {
      1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1588, 1784, 2000, 2238, 2500, 3500, 3832, 4252, 4713, 5220, 5777, 6389, 7000
   };

   const  int NPFBINS = 200;
   double pfrange[NPFBINS+1];
   for (int i = 0; i < NPFBINS+1; ++i) {
      pfrange[ i ] = i;
   }

   const int NPROBBINS = 100;
   double probrange[NPROBBINS+1];
   for ( int i = 0; i < NPROBBINS+1; ++i ) {
      probrange[i] = i/100;
   }

   #ifdef RECREATE_WEIGHTS
   FillWeightHistos(NBINS, ptrange, NPFBINS, pfrange);
   #endif

   TFile* f = TFile::Open("weightHistos.root");
   TH2D* gHist2;
   TH2D* qHist2;
   TH2D* gHist2_smooth;
   TH2D* qHist2_smooth;
   f->GetObject("gluon_nGenPF_probs", gHist2);
   f->GetObject("quark_nGenPF_probs", qHist2);
   f->GetObject("gluon_probs_smooth", gHist2_smooth);
   f->GetObject("quark_probs_smooth", qHist2_smooth);

   #ifdef SINGLE_TREE
   TFile* file = new TFile("outputGluonHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("outputGluonHistos.root", "recreate");
   #endif

   // Gluons
   TProfile* gProf1 = new TProfile("gluon_pt_resp", "", NBINS, ptrange);
   TProfile* gProf2 = new TProfile("gluon_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* gProf3 = new TProfile("gluon_nGenPF_prof", "", NBINS, ptrange);
   TProfile* gProf4 = new TProfile("gluon_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* gProf5 = new TProfile("gluon_pt_resp_nGenJetPF_w", "", NBINS, ptrange);

   TH2D* gHist3 = new TH2D("gluon_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* gHist4 = new TH1D("gluons_per_pt_bin", "", NBINS, ptrange);
   TH2D* gHist5 = new TH2D("gluon_genJetMass", "", NBINS, ptrange, NPFBINS, pfrange);
   TH2D* gHist6 = new TH2D("gluon_jetGirth", "", NBINS, ptrange, 100, 0., 0.5);
//   TH2D* gHist6 = new TH2D("gluon_")

   // Quarks
   TProfile* qProf1 = new TProfile("quark_pt_resp", "", NBINS, ptrange);
   TProfile* qProf2 = new TProfile("quark_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* qProf3 = new TProfile("quark_nGenPF_prof", "", NBINS, ptrange);
   TProfile* qProf4 = new TProfile("quark_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* qProf5 = new TProfile("quark_pt_resp_nGenJetPF_w", "", NBINS, ptrange);

   TH2D* qHist3 = new TH2D("quark_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* qHist4 = new TH1D("quarks_per_pt_bin", "", NBINS, ptrange);
   TH2D* qHist5 = new TH2D("quark_genJetMass", "", NBINS, ptrange, NPFBINS, pfrange);
   TH2D* qHist6 = new TH2D("quark_jetGirth", "", NBINS, ptrange, 100, 0., 0.5);

   //From weightHistos.root to outputHIstos.root
   TH2D* gHist2_clone = (TH2D*)gHist2->Clone("gluon_nGenPF_probs");
   TH2D* qHist2_clone = (TH2D*)qHist2->Clone("quark_nGenPF_probs");
   TH2D* gHist2_smooth_clone = (TH2D*)gHist2_smooth->Clone("gluon_nGenPF_probs_smooth");
   TH2D* qHist2_smooth_clone = (TH2D*)qHist2_smooth->Clone("quark_nGenPF_probs_smooth");

   double w;
   double gprob;
   double qprob;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
         gprob = gHist2->GetBinContent(gHist2_smooth->FindBin(jetPt, nGenJetPF));
         qprob = qHist2->GetBinContent(qHist2_smooth->FindBin(jetPt, nGenJetPF));
            if (isPhysG) {
               gProf1->Fill(jetPt, jetPt/genJetPt);
               gProf2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
               gProf3->Fill(jetPt, nGenJetPF);
               gHist4->Fill(jetPt);
               gHist5->Fill(jetPt, genJetMass);
               gHist6->Fill(jetPt, jetGirth);
               if (gprob > 0 && qprob > 0) {
                  w = qprob/gprob;
                  gProf5->Fill(jetPt, jetPt/genJetPt, w);
               }
               else w = 0;
               gHist3->Fill(jetPt, nGenJetPF, w);
               gProf4->Fill(jetPt, nGenJetPF, w);
            }
            else {
               if (isPhysUDS) {
                  qProf1->Fill(jetPt, jetPt/genJetPt);
                  qProf2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
                  qProf3->Fill(jetPt, nGenJetPF);
                  qHist4->Fill(jetPt);
                  qHist5->Fill(jetPt, genJetMass);
                  qHist6->Fill(jetPt, jetGirth);
                  if (gprob > 0 && qprob > 0) {
                     w = gprob/qprob;
                     qProf5->Fill(jetPt, jetPt/genJetPt, w);
                  }
                  else w = 0;
                  qHist3->Fill(jetPt, nGenJetPF, w);
                  qProf4->Fill(jetPt, nGenJetPF, w);
               }
            }
   }
}

file->Write();
file->Close();

f->Close();

}


void GluonHistosFill::FillWeightHistos(int nptbins, double* ptrange, int npfbins, double* pfrange)
{

   TFile* file = new TFile("./weightHistos.root", "recreate");

   TH2D* gHist1 = new TH2D("gluon_nGenPF_hist", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* gHist2 = new TH2D("gluon_nGenPF_probs", "", nptbins, ptrange, npfbins, pfrange);

   TH2D* qHist1 = new TH2D("quark_nGenPF_hist", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* qHist2 = new TH2D("quark_nGenPF_probs", "", nptbins, ptrange, npfbins, pfrange);   

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
            if (isPhysG) {
               gHist1->Fill(jetPt, nGenJetPF);
            }
            else {
               if (isPhysUDS) {
                  qHist1->Fill(jetPt, nGenJetPF);
               }
            }
      }
   }

   int gPFs;
   int qPFs;

   for (int xb = 1; xb != gHist1->GetNbinsX()+1; ++xb) {
      gPFs = gHist1->Integral(xb, xb, 0, 100);
      qPFs = qHist1->Integral(xb, xb, 0, 100);
      if (gPFs > 0 || qPFs > 0) { 
         for (int yb = 1; yb != gHist1->GetNbinsY()+1; ++yb) {
            if (gPFs > 0) {
               gHist2->SetBinContent(xb, yb, gHist1->GetBinContent(xb, yb) / gPFs);
            }
            if (qPFs > 0) {
               qHist2->SetBinContent(xb, yb, qHist1->GetBinContent(xb, yb) / qPFs);
            }
         }
      }
   }

   TH2D* gprobs_smooth = (TH2D*)gHist2->Clone("gluon_probs_smooth");
//   gprobs_smooth->SetAxisRange(2000, 3500, "X");
   gprobs_smooth->Smooth(1, "k5b");

   TH2D* qprobs_smooth = (TH2D*)qHist2->Clone("quark_probs_smooth");
//   qprobs_smooth->SetAxisRange(2000, 3500, "X");
   qprobs_smooth->Smooth(1, "k5b");

   cout << "probability calculated" << endl;

   TFile* f = TFile::Open("outputGluonHistos.root");

   file->Write();
   file->Close();
}