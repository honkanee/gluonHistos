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

   const int NMASSBINS = 200;
   double massrange[NMASSBINS+1];
   for (int i = 0; i < NMASSBINS+1; ++i) {
      massrange[ i ] = i;
   }

   #ifdef RECREATE_WEIGHTS
   FillWeightHistos(NBINS, ptrange, NPFBINS, pfrange, NMASSBINS, massrange);
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
   TProfile* gPtResp = new TProfile("gluon_pt_resp", "", NBINS, ptrange);
   TProfile* gPtRespNGenPF = new TProfile("gluon_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* gNGenPF = new TProfile("gluon_nGenPF_prof", "", NBINS, ptrange);
   TProfile* gNGenPF_w = new TProfile("gluon_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* gPtRespNGenPF_w = new TProfile("gluon_pt_resp_nGenJetPF_w", "", NBINS, ptrange);

   TH2D* gNGenPF_w_hist = new TH2D("gluon_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* gPerPtBin = new TH1D("gluons_per_pt_bin", "", NBINS, ptrange);

   // Quarks
   TProfile* qPtResp = new TProfile("quark_pt_resp", "", NBINS, ptrange);
   TProfile* qPtRespNGenPF = new TProfile("quark_pt_resp_nGenJetPF", "", NBINS, ptrange);
   TProfile* qNGenPF = new TProfile("quark_nGenPF_prof", "", NBINS, ptrange);
   TProfile* qNGenPF_w = new TProfile("quark_nGenPF_prof_w", "", NBINS, ptrange);
   TProfile* qPtRespNGenPF_w = new TProfile("quark_pt_resp_nGenJetPF_w", "", NBINS, ptrange);

   TH2D* qNGenPF_w_hist = new TH2D("quark_nGenPF_hist_w", "", NBINS, ptrange, NPFBINS, pfrange);
   TH1D* qPerPtBin = new TH1D("quarks_per_pt_bin", "", NBINS, ptrange);

   TH2D* nGenJetPF_probsRatio = new TH2D("nGenJetPF_likelyhood", "", NBINS, ptrange, NPFBINS, pfrange);

   TH2D* nGenPF_probs_sum = (TH2D*)gHist2_smooth->Clone("nGenPF_probs_sum");

//   nGenPF_probs_sum->Add(qHist2_smooth);
   nGenJetPF_probsRatio->Divide(qHist2_smooth, nGenPF_probs_sum);

   for (int xb = 1; xb != nGenJetPF_probsRatio->GetNbinsX()+1; ++xb) {
      for(int yb = 1; yb != nGenJetPF_probsRatio->GetNbinsY()+1; ++yb) {
            //nGenJetPF_probsRatio->SetBinContent(xb, yb, (nGenJetPF_probsRatio->GetBinContent(xb, yb)-0.5)*2);
            nGenJetPF_probsRatio->SetBinContent(xb, yb, nGenJetPF_probsRatio->GetBinContent(xb, yb)-1);

      }
   }

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
               gPtResp->Fill(jetPt, jetPt/genJetPt);
               gPtRespNGenPF->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
               gNGenPF->Fill(jetPt, nGenJetPF);
               gPerPtBin->Fill(jetPt);
               if (gprob > 0 && qprob > 0) {
                  w = qprob/gprob;
                  gPtRespNGenPF_w->Fill(jetPt, jetPt/genJetPt, w);
               }
               else w = 0;
               gNGenPF_w_hist->Fill(jetPt, nGenJetPF, w);
               gNGenPF_w->Fill(jetPt, nGenJetPF, w);
            }
            else {
               if (isPhysUDS) {
                  qPtResp->Fill(jetPt, jetPt/genJetPt);
                  qPtRespNGenPF->Fill(jetPt, jetPt/genJetPt, 1./nGenJetPF);
                  qNGenPF->Fill(jetPt, nGenJetPF);
                  qPerPtBin->Fill(jetPt);
                  if (gprob > 0 && qprob > 0) {
                     w = gprob/qprob;
                     qPtRespNGenPF_w->Fill(jetPt, jetPt/genJetPt, w);
                  }
                  else w = 0;
                  qNGenPF_w_hist->Fill(jetPt, nGenJetPF, w);
                  qNGenPF_w->Fill(jetPt, nGenJetPF, w);
               }
            }
   }
}

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

   //Quarks
   TH2D* qNGenPF = new TH2D("quark_nGenPF_hist", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* qGenJetMass = new TH2D("quark_genJetMass", "", nptbins, ptrange, nMassBins, pfrange);
   TH2D* qJetGirth = new TH2D("quark_jetGirth", "", nptbins, ptrange, 100, 0., 0.5);

   TH2D* qNGenPF_probs = new TH2D("quark_nGenPF_probs", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* qGenJetMass_probs = new TH2D("quark_genJetMass_probs", "", nptbins, ptrange, nMassBins, massrange);
   TH2D* qJetGirth_probs = new TH2D("quarks_jetGirth_probs", "", nptbins, ptrange, 100, 0., 0.5);

   int nGenPFsum;
   double nUEPF_per_A;
   double nUEPF;
   double dR_in = 0.8;
   double dR_out = 1.0;
   TH2D* UEPFs = new TH2D("UE_nPF", "", nptbins, ptrange, npfbins, pfrange);
   TH2D* UEPF_probs = new TH2D("UE_nPF_probs", "", nptbins, ptrange, npfbins, pfrange);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {
         if (isPhysG) {
            gNGenPF->Fill(jetPt, nGenJetPF);
            gGenJetMass->Fill(jetPt, genJetMass);
            gJetGirth->Fill(jetPt, jetGirth);
            nGenPFsum = 0;
            for (Float_t p : PF_dR) {
               if (p > dR_in && p < dR_out) {
                  ++nGenPFsum;
               }
            }
            nUEPF_per_A = nGenPFsum / (M_PI * (dR_out*dR_out - dR_in*dR_in));
            nUEPF = jetArea*nUEPF_per_A;
            UEPFs->Fill(jetPt, nGenPFsum);
         }
         else {
            if (isPhysUDS) {
               qNGenPF->Fill(jetPt, nGenJetPF);
               qGenJetMass->Fill(jetPt, genJetMass);
               qJetGirth->Fill(jetPt, jetGirth);
               nGenPFsum = 0;
               for (auto p : PF_dR) {
                  if (p > dR_in && p < dR_out) {
                     ++nGenPFsum;
                  }
               }
               nUEPF_per_A = nGenPFsum / (M_PI * (dR_out*dR_out - dR_in*dR_in));
               nUEPF = jetArea*nUEPF_per_A;
               UEPFs->Fill(jetPt, nGenPFsum);
            }
         }
      }
   }

   CalculateProbs(UEPFs, nullptr, UEPF_probs, nullptr);

   gNGenPF->Smooth(1, "k5b");
   qNGenPF->Smooth(1, "k5b");

   CalculateProbs(gNGenPF, qNGenPF, gNGenPF_probs, qNGenPF_probs);
   CalculateProbs(gGenJetMass, qGenJetMass, gGenJetMass_probs, qGenJetMass_probs);
   CalculateProbs(gJetGirth, qJetGirth, gJetGirth_probs, qJetGirth_probs);

   TH2D* gprobs_smooth = (TH2D*)gNGenPF_probs->Clone("gluon_probs_smooth");
//   gprobs_smooth->SetAxisRange(2000, 3500, "X");
//   gprobs_smooth->Smooth(1, "k5b");

   TH2D* qprobs_smooth = (TH2D*)qNGenPF_probs->Clone("quark_probs_smooth");
//   qprobs_smooth->SetAxisRange(2000, 3500, "X");
//   qprobs_smooth->Smooth(1, "k5b");

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