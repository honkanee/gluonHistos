#define GluonHistosFill_cxx

#include "GluonHistosFill.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

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

   #ifdef SINGLE_TREE
   TFile* file = new TFile("outputGluonHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("outputGluonHistos.root", "recreate");   
   #endif

   const int NBINS = 79;
   const double ptrange[NBINS+1] = {
      1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
   };

   // Gluons
   TProfile* gHist1 = new TProfile("gluon_pt_resp", "", NBINS, ptrange);
   TProfile* gHist2 = new TProfile("gluon_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* gHist3 = new TProfile("gluon_nGenPF", "", NBINS, ptrange);
   TH1D* gHist4 = new TH1D("gluon_nGenPF_w", "", NBINS, ptrange);

   // Quarks
   TProfile* qHist1 = new TProfile("quark_pt_resp", "", NBINS, ptrange);
   TProfile* qHist2 = new TProfile("quark_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* qHist3 = new TProfile("quark_nGenPF", "", NBINS, ptrange);
   TH1D* qHist4 = new TH1D("quark_nGenPF_w", "", NBINS, ptrange);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
            if (isPhysG) {
               gHist1->Fill(jetPt, jetPt/genJetPt);
               gHist2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
               gHist3->Fill(jetPt, nGenJetPF);
            }
            else {
               qHist1->Fill(jetPt, jetPt/ genJetPt);
               qHist2->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
               qHist3->Fill(jetPt, nGenJetPF);
            }
      };
   }

   TH1D* q_g_nGenPF_ratio = new TH1D("q_g_nGenPF_ratio", "", NBINS, ptrange);
   q_g_nGenPF_ratio->Divide(qHist3->ProjectionX(""), gHist3->ProjectionX(""));

   gHist4->Multiply(gHist3->ProjectionX(""), q_g_nGenPF_ratio);
   qHist4->Divide(qHist3->ProjectionX(""), q_g_nGenPF_ratio);

   file->Write();
   file->Close();
}
