#define GluonHistosFill_cxx

#include "GluonHistosFill.h"

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

#include <vector>
#include <ctime>
#include <chrono>

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

auto start = chrono::system_clock::now(); //timer

   if (fChain == 0) return;

   #ifdef SINGLE_TREE
   TFile* file = new TFile("outputGluonHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("outputGluonHistos.root", "recreate");   
   #endif

   const int NBINS = 79;
   const double ptrange[NBINS+1] = {
      1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
   };

   const int N2BINS = 100;
   double pfrange[N2BINS+1];
   double probrange[N2BINS+1];     
   for ( int i = 0; i < N2BINS+1; ++i ) {
      pfrange[ i ] = i;
      probrange[i] = i/100;
   }

   // Gluons
   TProfile* gProf1 = new TProfile("gluon_pt_resp", "", NBINS, ptrange);
   TProfile* gProf2 = new TProfile("gluon_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* gProf3 = new TProfile("gluon_nGenPF_prof", "", NBINS, ptrange);
   TProfile* gProf4 = new TProfile("gluon_nGenPF_prof_w", "", NBINS, ptrange);

   TH2D* gHist1 = new TH2D("gluon_nGenPF_hist", "", NBINS, ptrange, N2BINS, pfrange);
   TH2D* gHist2 = new TH2D("gluon_nGenPF_probs", "", NBINS, ptrange, N2BINS, pfrange);
   TH2D* gHist3 = new TH2D("gluon_nGenPF_hist_w", "", NBINS, ptrange, N2BINS, pfrange);

   // Quarks
   TProfile* qProf1 = new TProfile("quark_pt_resp", "", NBINS, ptrange);
   TProfile* qProf2 = new TProfile("quark_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* qProf3 = new TProfile("quark_nGenPF_prof", "", NBINS, ptrange);
   TProfile* qProf4 = new TProfile("quark_nGenPF_prof_w", "", NBINS, ptrange);

   TH2D* qHist1 = new TH2D("quark_nGenPF_hist", "", NBINS, ptrange, N2BINS, pfrange);
   TH2D* qHist2 = new TH2D("quark_nGenPF_probs", "", NBINS, ptrange, N2BINS, pfrange);
   TH2D* qHist3 = new TH2D("quark_nGenPF_hist_w", "", NBINS, ptrange, N2BINS, pfrange);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
            if (isPhysG) {
               gProf1->Fill(jetPt, jetPt/genJetPt);
               gProf2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
               gProf3->Fill(jetPt, nGenJetPF);
               gHist1->Fill(jetPt, nGenJetPF);
            }
            else {
               if (isPhysUDS) {
                  qProf1->Fill(jetPt, jetPt/ genJetPt);
                  qProf2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
                  qProf3->Fill(jetPt, nGenJetPF);
                  qHist1->Fill(jetPt, nGenJetPF);
               }
            }
      };
   }
   cout << "first loop done" << endl;

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

   cout << "probability calculated" << endl;

   double w;
   double gprob;
   double qprob;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
         w = 1;
         gprob = gHist2->GetBinContent(jetPt, nGenJetPF);
         qprob = qHist2->GetBinContent(jetPt, nGenJetPF);
            if (isPhysG) {
               if (gprob > 0) {
                  w = qprob/gprob;
               }
               gHist3->Fill(jetPt, nGenJetPF, w);
               gProf4->Fill(jetPt, nGenJetPF, w);
            }
            else {
               if (isPhysUDS) {
                  if (qprob > 0) {
                     w = gprob/qprob;
                  }
                  qHist3->Fill(jetPt, nGenJetPF, w);
                  qProf4->Fill(jetPt, nGenJetPF, w);
               }
            }
   };
}

file->Write();
file->Close();


// timer
auto end = std::chrono::system_clock::now();
std::chrono::duration<double> elapsed_sec = end-start;
std::time_t end_time = std::chrono::system_clock::to_time_t(end);
std::cout << "finished at " << std::ctime(&end_time)
         << "elapsed time: " << elapsed_sec.count()/60 << "min";
}
