#define GluonHistosFill_cxx

#include "GluonHistosFill.h"

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TF1.h>

#include <vector>

char* char_add(string a, string b) {
   string str = a + b;
   char *cstr = new char[str.length() + 1];
   strcpy(cstr, str.c_str());
   return cstr;
}

enum {GLUON, QUARK, ALL};

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
fChain->SetBranchStatus("jetRawPt", 1);

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

   TH2D* gJetGirth_probs;
   TH2D* qJetGirth_probs;

   f->GetObject("gluon_nGenPF_probs", gHist2);
   f->GetObject("quark_nGenPF_probs", qHist2);
   f->GetObject("UE_nPF_perA", nUEPF_perA_hist);
   f->GetObject("gluon_nGenPFSig_probs", gHist2_sig);
   f->GetObject("quark_nGenPFSig_probs", qHist2_sig);
   f->GetObject("gluon_jetGirth_probs", gJetGirth_probs);
   f->GetObject("quark_jetGirth_probs", qJetGirth_probs);

   #ifdef SINGLE_TREE
   TFile* file = new TFile("outputGluonHistos_single.root", "recreate");
   #else
   TFile* file = new TFile("outputGluonHistos.root", "recreate");
   #endif

   vector<TProfile*> ptResp;
   vector<TProfile*> nGenPF;
   vector<TProfile*> nGenPF_w;
   vector<TProfile*> ptRespNGenPF_w;
   vector<TProfile*> ptRespNGenPFSig_w;

   vector<TProfile*> ptRespJetGirth_w;

   vector<TH2D*> nGenPF_w_hist;
   vector<TH1D*> perPtBin;
   vector<TH2D*> ptResp_hist;
   vector<TH2D*> ptResp_hist_norm;
   vector<TH2D*> ptRespNGenPF_w_hist;
   vector<TH2D*> ptRespNGenPFSig_w_hist;

   vector<TH1D*> ptResp_gaus;
   vector<TH1D*> ptRespNGenPF_w_gaus;
   vector<TH1D*> ptRespNGenPFSig_w_gaus;

   for (string str : {"gluon_", "quark_"}) {
      ptResp.push_back(new TProfile(char_add(str, "pt_resp"), "", NBINS, ptrange));
      nGenPF.push_back(new TProfile(char_add(str, "nGenPF_prof"), "", NBINS, ptrange));
      nGenPF_w.push_back(new TProfile(char_add(str, "nGenPF_prof_w"), "", NBINS, ptrange));
      ptRespNGenPF_w.push_back(new TProfile(char_add(str, "pt_resp_nGenJetPF_w"), "", NBINS, ptrange));
      ptRespNGenPFSig_w.push_back(new TProfile(char_add(str, "pt_resp_nGenJetPFSig_w"), "", NBINS, ptrange));

      ptRespJetGirth_w.push_back(new TProfile(char_add(str, "pt_resp_jetGirth_w"), "", NBINS, ptrange));

      nGenPF_w_hist.push_back(new TH2D(char_add(str, "nGenPF_hist_w"), "", NBINS, ptrange, NPFBINS, pfrange));
      perPtBin.push_back(new TH1D(char_add(str, "gluons_per_pt_bin"), "", NBINS, ptrange));
      ptResp_hist.push_back(new TH2D(char_add(str, "pt_resp_hist"), "", NBINS, ptrange, 100, 0.5, 1.5));
      ptResp_hist_norm.push_back(new TH2D(char_add(str, "pt_resp_hist_norm"), "", NBINS, ptrange, 100, 0.5,1.5));
      ptRespNGenPF_w_hist.push_back(new TH2D(char_add(str, "pt_resp_nGenJetPF_w_hist"), "", NBINS, ptrange, 100, 0.5, 1.5));
      ptRespNGenPFSig_w_hist.push_back(new TH2D(char_add(str, "pt_resp_nGenJetPFSig_w_hist"), "", NBINS, ptrange, 100, 0.5, 1.5));

      ptResp_gaus.push_back(new TH1D(char_add(str, "pt_resp_gaus"), "", NBINS, ptrange));
      ptRespNGenPF_w_gaus.push_back(new TH1D(char_add(str, "pt_resp_nGenJetPF_w_gaus"), "", NBINS, ptrange));
      ptRespNGenPFSig_w_gaus.push_back(new TH1D(char_add(str, "pt_resp_nGenJetPFSig_w_gaus"), "", NBINS, ptrange));
   }

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
   TH2D* PtResp_hist = new TH2D("pt_resp_hist", "", NBINS, ptrange, 100, 0.5, 1.5);
   TH2D* PtRespNGenPF_wg_hist = new TH2D("pt_resp_nGenJetPF_wg_hist", "", NBINS, ptrange, 100, 0.5, 1.5);
   TH2D* PtRespNGenPFSig_wg_hist = new TH2D("pt_resp_nGenJetPFSig_wg_hist", "", NBINS, ptrange, 100, 0.5, 1.5);
   TH2D* PtRespNGenPF_wq_hist = new TH2D("pt_resp_nGenJetPF_wq_hist", "", NBINS, ptrange, 100, 0.5, 1.5);
   TH2D* PtRespNGenPFSig_wq_hist = new TH2D("pt_resp_nGenJetPFSig_wq_hist", "", NBINS, ptrange, 100, 0.5, 1.5);

   TH1D* PtResp_gaus = new TH1D("pt_resp_gaus", "", NBINS, ptrange);
   TH1D* PtRespNGenPF_wg_gaus = new TH1D("pt_resp_nGenJetPF_wg_gaus", "", NBINS, ptrange);
   TH1D* PtRespNGenPFSig_wg_gaus = new TH1D("pt_resp_nGenJetPFSig_wg_gaus", "", NBINS, ptrange);
   TH1D* PtRespNGenPF_wq_gaus = new TH1D("pt_resp_nGenJetPF_wq_gaus", "", NBINS, ptrange);
   TH1D* PtRespNGenPFSig_wq_gaus = new TH1D("pt_resp_nGenJetPFSig_wq_gaus", "", NBINS, ptrange);

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
   double nGenJetPF_w;
   double jetGirth_w;
   double gprob;
   double qprob;
   double gprob_jetGirth;
   double qprob_jetGirth;

   double nUEPF_per_A = nUEPF_perA_hist->GetMean();
   int nPFsum;
   double dR_in = 0.8;
   double dR_out = 1.0;
   double nUEPF;
   double nGenPFSig;

   double w_sig;
   double gprob_sig;
   double qprob_sig;

   int type;

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
            if (rawWeighted) {
               usedWeightPt = jetRawPt;
            }
            else {
               usedWeightPt = genJetPt;
            }
         }

         resp = jetPt/genJetPt; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

         gprob = gHist2->GetBinContent(gHist2->FindBin(usedWeightPt, nGenJetPF));
         qprob = qHist2->GetBinContent(qHist2->FindBin(usedWeightPt, nGenJetPF));

         gprob_jetGirth = gJetGirth_probs->GetBinContent(gJetGirth_probs->FindBin(usedWeightPt, jetGirth));
         qprob_jetGirth = qJetGirth_probs->GetBinContent(qJetGirth_probs->FindBin(usedWeightPt, jetGirth));

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
            type = GLUON;
            if (gprob_sig > 0 && qprob_sig > 0) {
               w_sig = qprob_sig/gprob_sig;
            }
            else w_sig = 0;
            if (gprob > 0 && qprob > 0) {
               nGenJetPF_w = qprob/gprob;
            }
            else nGenJetPF_w = 0;
            if (gprob_jetGirth > 0 && qprob_jetGirth > 0) {
               jetGirth_w = qprob_jetGirth/gprob_jetGirth;
            }
            else jetGirth_w = 0;
         }
         else {
            if (isPhysUDS) {
               type = QUARK;
               if (gprob_sig > 0 && qprob_sig > 0) {
                  w_sig = gprob_sig/qprob_sig;
               }
               else w_sig = 0;
               if (gprob > 0 && qprob > 0) {
                  nGenJetPF_w = gprob/qprob;
               }
               else nGenJetPF_w = 0;
               if (gprob_jetGirth > 0 && qprob_jetGirth > 0) {
                  jetGirth_w = gprob_jetGirth/qprob_jetGirth;
               }
               else jetGirth_w = 0;
            }
         }
         if (type == GLUON || type == QUARK) {
            type = QUARK;
            ptResp.at(type)->Fill(usedPlotPt, resp);
            ptResp_hist.at(type)->Fill(usedPlotPt, resp);
            nGenPF.at(type)->Fill(usedPlotPt, nGenJetPF);
            perPtBin.at(type)->Fill(usedPlotPt);
            ptRespNGenPFSig_w.at(type)->Fill(usedPlotPt, resp, w_sig);
            ptRespNGenPFSig_w_hist.at(type)->Fill(usedPlotPt, resp, w_sig);
            ptRespNGenPF_w.at(type)->Fill(usedPlotPt, resp, nGenJetPF_w);
            ptRespNGenPF_w_hist.at(type)->Fill(usedPlotPt, resp, nGenJetPF_w);
            ptRespJetGirth_w.at(type)->Fill(usedPlotPt, resp, jetGirth_w);
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

   double fit_up = 1.2;
   double fit_down = 0.8;

   double cap = 0.12;

   double mean;
   TH1D* py;
   TF1 *g1;
   TH2D* gbounds = new TH2D("gluon_resp_gaus_bounds", "", NBINS, ptrange, 100, 0.5, 1.5);
   TH2D* qbounds = new TH2D("quark_resp_gaus_bounds", "", NBINS, ptrange, 100, 0.5, 1.5);

   vector<TH2D*> respHists = {PtResp_hist, PtRespNGenPF_wg_hist, PtRespNGenPFSig_wg_hist, PtRespNGenPFSig_wq_hist, PtRespNGenPFSig_wq_hist,
                              ptResp_hist.at(GLUON), ptRespNGenPF_w_hist.at(GLUON), ptRespNGenPFSig_w_hist.at(GLUON),
                              ptResp_hist.at(QUARK), ptRespNGenPF_w_hist.at(QUARK), ptRespNGenPFSig_w_hist.at(QUARK)};
   vector<TH1D*> gausHists = {PtResp_gaus, PtRespNGenPF_wg_gaus, PtRespNGenPFSig_wg_gaus, PtRespNGenPFSig_wq_gaus, PtRespNGenPFSig_wq_gaus,
                              ptResp_gaus.at(GLUON), ptRespNGenPF_w_gaus.at(GLUON), ptRespNGenPFSig_w_gaus.at(GLUON),
                              ptResp_gaus.at(QUARK), ptRespNGenPF_w_gaus.at(QUARK), ptRespNGenPFSig_w_gaus.at(QUARK)};

   for (int i = 0; i != respHists.size(); ++i) {
      for (int xb = 1; xb != respHists.at(i)->GetNbinsX()+1; ++xb) {
         respHists.at(i)->GetXaxis()->SetRange(xb, xb);
         mean = respHists.at(i)->GetMean(2);
         g1 = new TF1("g1","gaus",mean - cap, mean + cap);
         py = respHists.at(i)->ProjectionY("_py", xb, xb);
         py->Fit(g1, "Q", "");
         gausHists.at(i)->SetBinContent(xb, g1->GetParameter("Mean"));
         
         if (respHists.at(i) == ptResp_hist.at(GLUON)) {
            gbounds->Fill(respHists.at(i)->GetXaxis()->GetBinLowEdge(xb), mean + cap);
            gbounds->Fill(respHists.at(i)->GetXaxis()->GetBinLowEdge(xb), mean - cap);
         }
         if (respHists.at(i) == ptResp_hist.at(QUARK)) {
            qbounds->Fill(respHists.at(i)->GetXaxis()->GetBinLowEdge(xb), mean + cap);
            qbounds->Fill(respHists.at(i)->GetXaxis()->GetBinLowEdge(xb), mean - cap);
         }
      }
   }

   CalculateProbs(ptResp_hist.at(GLUON), ptResp_hist.at(QUARK), ptResp_hist_norm.at(GLUON), ptResp_hist_norm.at(QUARK));

   const int nbins[3] = {5, 5, 5};
   const double mins[3] = {0, 0, 0};
   const double maxs[3] = {5, 5, 5};

   THnD* testi = new THnD("testi", "", 3, nbins, mins, maxs);
   THnD* testi_probs = new THnD("testi_probs", "", 3, nbins, mins, maxs);

   double x[3] = {0, 0, 0};
   testi->Fill(x);

   int y[3] = {1, 1, 1};
   cout << "a: " << testi->GetBinContent(y) << endl;

   CalculateProbsNdim(testi, testi_probs);


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

   vector<TH2D*> nGenPF;
   vector<TH2D*> genJetMass_histos;
   vector<TH2D*> jetGirth_histos;

   vector<TH2D*> nGenPF_probs;
   vector<TH2D*> genJetMass_probs;
   vector<TH2D*> jetGirth_probs;
   vector<TH2D*> nGenPFSig;
   vector<TH2D*> nGenPFSig_probs;
   vector<TProfile*> nGenPFSig_prof;

   vector<TH1D*> dR_nPF_fromPV0;
   vector<TH1D*> dR_nPF_fromPV1;
   vector<TH1D*> dR_nPF_fromPV2;
   vector<TH1D*> dR_nPF_fromPV3;

   for (string str : {"gluon_", "quark_"}) {
      nGenPF.push_back(new TH2D(char_add(str, "nGenPF_hist"), "", nptbins, ptrange, npfbins, pfrange));
      genJetMass_histos.push_back(new TH2D(char_add(str, "genJetMass"), "", nptbins, ptrange, nMassBins, massrange));
      jetGirth_histos.push_back(new TH2D(char_add(str, "jetGirth"), "", nptbins, ptrange, 100, 0., 0.5));

      nGenPF_probs.push_back(new TH2D(char_add(str, "nGenPF_probs"), "", nptbins, ptrange, npfbins, pfrange));
      genJetMass_probs.push_back(new TH2D(char_add(str, "genJetMass_probs"), "", nptbins, ptrange, nMassBins, massrange));
      jetGirth_probs.push_back(new TH2D(char_add(str, "jetGirth_probs"), "", nptbins, ptrange, 100, 0., 0.5));
      nGenPFSig.push_back(new TH2D(char_add(str, "nGenPFSig"), "", nptbins, ptrange, 300, -50, 250));
      nGenPFSig_probs.push_back(new TH2D(char_add(str, "nGenPFSig_probs"), "", nptbins, ptrange, 300, -50, 250));
      nGenPFSig_prof.push_back(new TProfile(char_add(str, "nGenPFSig_prof"), "", nptbins, ptrange));

      dR_nPF_fromPV0.push_back(new TH1D(char_add(str, "dR_nPF_fromPV0"), "", 100, 0, 1.5));
      dR_nPF_fromPV1.push_back(new TH1D(char_add(str, "dR_nPF_fromPV1"), "", 100, 0, 1.5));
      dR_nPF_fromPV2.push_back(new TH1D(char_add(str, "dR_nPF_fromPV2"), "", 100, 0, 1.5));
      dR_nPF_fromPV3.push_back(new TH1D(char_add(str, "dR_nPF_fromPV3"), "", 100, 0, 1.5));
   }

   int nPFsum;
   double nUEPF_per_A;
   double nUEPF;
   double dR_in = 0.8;
   double dR_out = 1.0;
   double nGenPFSig_value;

   TH1D* nUEPF_perA_hist = new TH1D("UE_nPF_perA", "", int(80./2.75)+3, 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_nGenPF_hist = new TH2D("nUEPF_perA_nGenPF_hist", "", npfbins, pfrange, int(80./2.75)+3, 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_nGenPF_probs = new TH2D("nUEPF_perA_nGenPF_probs", "", npfbins, pfrange, int(80./2.75)+3, 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_genJetPt_hist = new TH2D("nUEPF_perA_genJetPt_hist", "", nptbins, ptrange, int(80./2.75)+3, 0, int(80./2.75)*2.75);
   TH2D* nUEPF_perA_genJetPt_probs = new TH2D("nUEPF_perA_genJetPt_probs", "", nptbins, ptrange, int(80./2.75)+3, 0, int(80./2.75)*2.75);
   TProfile* nUEPF_nGenPF = new TProfile("nUEPF_nGenPF", "", npfbins, pfrange);

   TH1D* dR_nPF_fromPV0_all = new TH1D("dR_nPF_fromPV0", "", 100, 0, 1.5);
   TH1D* dR_nPF_fromPV1_all = new TH1D("dR_nPF_fromPV1", "", 100, 0, 1.5);
   TH1D* dR_nPF_fromPV2_all = new TH1D("dR_nPF_fromPV2", "", 100, 0, 1.5);
   TH1D* dR_nPF_fromPV3_all = new TH1D("dR_nPF_fromPV3", "", 100, 0, 1.5);

   int all_events = 0;
   int gluon_jets = 0;
   int quark_jets = 0;
   vector<int> PV0_jets = {0, 0, 0};
   vector<int> PV1_jets = {0, 0, 0};
   vector<int> PV2_jets = {0, 0, 0};
   vector<int> PV3_jets = {0, 0, 0};

   int type;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (fabs(jetEta)<1.3 && jetPtOrder<2 && fpclassify(genJetPt) == FP_NORMAL && genJetPt > 0) {

         ++all_events;
         if (isPhysG) {
            type = GLUON;
            ++gluon_jets;
         }
         else {
            if (isPhysUDS) {
               type = QUARK;
               ++quark_jets;
            }
         }
         if (recoWeighted) {
            usedWeightPt = jetPt;
         }
         else {
            if (rawWeighted) {
               usedWeightPt = jetRawPt;
            } 
            else {
               usedWeightPt = genJetPt;
            }
         }
         for (int i = 0; i != nPF; ++i) {
            if (PF_fromPV[i] == 0) {
               ++PV0_jets.at(ALL);
               dR_nPF_fromPV0_all->Fill(PF_dR[i]);
               if (isPhysG || isPhysUDS) {
                  ++PV0_jets.at(type);
                  dR_nPF_fromPV0.at(type)->Fill(PF_dR[i]);
               }
            }
            if (PF_fromPV[i] == 1) {
               ++PV1_jets.at(ALL);
               dR_nPF_fromPV1_all->Fill(PF_dR[i]);
               if (isPhysG|| isPhysUDS) {
                  ++PV1_jets.at(type);
                  dR_nPF_fromPV1.at(type)->Fill(PF_dR[i]);
               }
            }
            if (PF_fromPV[i] == 2) {
               ++PV2_jets.at(ALL);
               dR_nPF_fromPV2_all->Fill(PF_dR[i]);
               if (isPhysG || isPhysUDS) {
                  ++PV2_jets.at(type);
                  dR_nPF_fromPV2.at(type)->Fill(PF_dR[i]);
               }
            }
            if (PF_fromPV[i] == 3) {
               ++PV3_jets.at(ALL);
               dR_nPF_fromPV3_all->Fill(PF_dR[i]);
               if (isPhysUDS || isPhysG) {
                  ++PV3_jets.at(type);
                  dR_nPF_fromPV3.at(type)->Fill(PF_dR[i]);
               }
            }
         }
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
         nUEPF = jetArea*nUEPF_per_A;
         nGenPFSig_value = nGenJetPF - nUEPF;

         if (type == GLUON || type == QUARK) {
            nGenPF.at(type)->Fill(usedWeightPt, nGenJetPF);
            genJetMass_histos.at(type)->Fill(usedWeightPt, genJetMass);
            jetGirth_histos.at(type)->Fill(usedWeightPt, jetGirth);
            if (isPFSig_perEvent) {
               nGenPFSig.at(type)->Fill(usedWeightPt, nGenPFSig_value);
            }
         }
      }
   }

   cout << "first loop done" << endl;

   double a;
   double low;
   double up;
   double normalized;

   vector<TH1D*> fromPV_histos = {dR_nPF_fromPV0_all, dR_nPF_fromPV1_all, dR_nPF_fromPV2_all, dR_nPF_fromPV3_all, 
                     dR_nPF_fromPV0.at(GLUON), dR_nPF_fromPV1.at(GLUON), dR_nPF_fromPV2.at(GLUON), dR_nPF_fromPV3.at(GLUON),
                     dR_nPF_fromPV0.at(QUARK), dR_nPF_fromPV1.at(QUARK), dR_nPF_fromPV2.at(QUARK), dR_nPF_fromPV3.at(QUARK)};
   vector<int> counts = {all_events, all_events, all_events, all_events, gluon_jets, gluon_jets, gluon_jets, gluon_jets, quark_jets, quark_jets, quark_jets, quark_jets};
   TH1D* hist;
   for (int i = 0; i != fromPV_histos.size(); ++i) {
      hist = fromPV_histos.at(i);
      for (int xb = 1; xb != hist->GetNbinsX()+1; ++xb) {
         low = hist->GetXaxis()->GetBinLowEdge(xb);
         up = hist->GetXaxis()->GetBinLowEdge(xb+1);
         a = M_PI * (up * up - low * low);
         normalized = hist->GetBinContent(xb) / a;
         hist->SetBinContent(xb, normalized);
      }
      hist->Scale(1./counts.at(i));
      cout << counts.at(i) << endl;
   }

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
            else {
               if (rawWeighted) {
                  usedWeightPt = jetRawPt;
               }
               else {
                  usedWeightPt = genJetPt;
               }
            } 

            if (isPhysG) {
               type = GLUON;
            }
            if (isPhysUDS) {
                  type = QUARK;
               }
            if (type == GLUON || type == QUARK) {
               nUEPF = jetArea*nUEPF_per_A;
               nGenPFSig_value = nGenJetPF - nUEPF;
               nGenPFSig.at(type)->Fill(usedWeightPt, nGenPFSig_value);
               nGenPFSig_prof.at(type)->Fill(usedWeightPt, nGenPFSig_value);
            }
         }
      }
   }

   nGenPF.at(GLUON)->Smooth(1, "k5b");
   nGenPFSig.at(GLUON)->Smooth(1, "k5b");
   nGenPF.at(QUARK)->Smooth(1, "k5b");
   nGenPFSig.at(QUARK)->Smooth(1, "k5b");

   CalculateProbs(nUEPF_perA_genJetPt_hist, nullptr, nUEPF_perA_genJetPt_probs, nullptr);
   CalculateProbs(nUEPF_perA_nGenPF_hist, nullptr, nUEPF_perA_nGenPF_probs, nullptr);

   CalculateProbs(nGenPF.at(GLUON), nGenPF.at(QUARK), nGenPF_probs.at(GLUON), nGenPF_probs.at(QUARK));
   CalculateProbs(nGenPFSig.at(GLUON), nGenPFSig.at(QUARK), nGenPFSig_probs.at(GLUON), nGenPFSig_probs.at(QUARK));

   CalculateProbs(genJetMass_histos.at(GLUON), genJetMass_histos.at(QUARK), genJetMass_probs.at(GLUON), genJetMass_probs.at(QUARK));
   CalculateProbs(jetGirth_histos.at(GLUON), jetGirth_histos.at(QUARK), jetGirth_probs.at(GLUON), jetGirth_probs.at(QUARK));

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

void GluonHistosFill::CalculateProbs3D(TH3D* gHist, TH3D* qHist, TH3D* gProbs, TH3D* qProbs)
{
   double gIntegral;
   double qIntegral;
   for (int xb = 1; xb != gHist->GetNbinsX()+1; ++xb) {
      gIntegral = gHist->Integral(xb, xb, 0, -1, 0, -1);
      if (qHist != nullptr) {
         qIntegral = qHist->Integral(xb, xb, 0, -1, 0, -1);
      }
      if (gIntegral > 0 || (qIntegral > 0 && qProbs != nullptr)) {
         for (int zb = 1; zb != gHist->GetNbinsZ()+1; ++zb) {
            for (int yb = 1; yb != gHist->GetNbinsY()+1; ++yb) {
               if (gIntegral > 0) {
                  gProbs->SetBinContent(xb, yb, zb, gHist->GetBinContent(xb, yb) / gIntegral);
               }
               if (qIntegral > 0 && qHist != nullptr) {
                  qProbs->SetBinContent(xb, yb, zb, qHist->GetBinContent(xb, yb) / qIntegral);
               }
            }
         }
      }
   }
}

void  GluonHistosFill::CalculateProbsNdim(THnD* hist, THnD* probs)
{
   int d;
   int D;
   double integral;
   vector<double> integrals;
   int dims = hist->GetNdimensions();
   int bin[dims];
   int i;
   bool start_over;

   for (int i = 0; i != dims; ++i) {
      bin[i] = 1;
   }

   for (int xb = 1; xb != hist->GetAxis(1)->GetNbins() +1; ++xb) {
      bin[0] = xb;
      for (int i = 1; i != dims; ++i) {
         bin[i] = 1;
      }
      integral = 0;
      D = 1;
      d =1;
      integral += hist->GetBinContent(bin);
      cout << bin[0] << bin[1] << bin[2] << endl;
      while(!(bin[dims-1] == hist->GetAxis(dims-1)->GetNbins() && bin[dims-2] == hist->GetAxis(dims-2)->GetNbins())) {
         while (bin[d] == hist->GetAxis(d)->GetNbins()) {
            ++d;
            for (int i = 1; i != d; ++i) {
               bin[i] = 1;
            }
            start_over = true;
         }
         if (start_over) {
            ++bin[d];
            d = 1;
            start_over = false;
         }
         else {
            ++bin[d];
         }
         integral += hist->GetBinContent(bin);
         cout << bin[0] << bin[1] << bin[2] << endl;
      }
      integrals.push_back(integral);
   }
      for (int xb = 1; xb != hist->GetAxis(1)->GetNbins() +1; ++xb) {
      bin[0] = xb;
      for (int i = 1; i != dims; ++i) {
         bin[i] = 1;
      }
      D = 1;
      d =1;
      integral = integrals.at(xb-1);
      if (integral != 0) {
         probs->SetBinContent(bin, hist->GetBinContent(bin)/integral);
         while(!(bin[dims-1] == hist->GetAxis(dims-1)->GetNbins() && bin[dims-2] == hist->GetAxis(dims-2)->GetNbins())) {
            while (bin[d] == hist->GetAxis(d)->GetNbins()) {
               ++d;
               for (int i = 1; i != d; ++i) {
                  bin[i] = 1;
               }
               start_over = true;
            }
            if (start_over) {
               ++bin[d];
               d = 1;
               start_over = false;
            }
            else {
               ++bin[d];
            }
            probs->SetBinContent(bin, hist->GetBinContent(bin)/integral);
         }
      }
   }
}