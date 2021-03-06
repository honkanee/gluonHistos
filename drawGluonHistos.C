#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <vector>

#include "tdrstyle_mod15.C"

using namespace std;

char* char_add(string a, string b) {
   string str = a + b;
   char *cstr = new char[str.length() + 1];
   strcpy(cstr, str.c_str());
   return cstr;
}

enum {GLUON, QUARK, ALL};
enum {G, UD, S};

void drawGluonHistos() {
    setTDRStyle();

    TFile* f = TFile::Open("outputGluonHistos.root");
    TFile* f2 = TFile::Open("weightHistos.root");

    //All
    TH2D* nUEPF_pt_hist;
    TH2D* nUEPF_nGenPF_hist;

    f2->GetObject("nUEPF_perA_genJetPt_probs", nUEPF_pt_hist);
    f2->GetObject("nUEPF_perA_nGenPF_probs", nUEPF_nGenPF_hist);

    vector<TProfile*> pt_resp(3);
    vector<TProfile*> pt_resp_nGenJetPF(3);
    vector<TProfile*> nGenPF(3);
    vector<TProfile*> nGenPF_w(3);
    vector<TProfile*> pt_resp_nGenJetPF_w(3);
    vector<TProfile*> pt_resp_nGenJetPFSig_w(3);
    vector<TProfile*> nGenPFSig_prof(3);
    vector<TProfile*> pt_resp_jetGirth_w(3);
    vector<TProfile*> pt_resp_genJetMass_w(3);
    vector<TProfile*> pt_resp_genJetLHA_w(3);
    vector<TProfile*> pt_resp_genJetWidth_w(3);
    vector<TProfile*> pt_resp_genJetThrust_w(3);
    vector<TProfile*> pt_resp_genJetMultiplicity_w(3);
    vector<TProfile*> pt_resp_genJetPtD_w(3);
    vector<TProfile*> pt_resp_nGenJetPF_genJetWidth_w(3);
    vector<TProfile*> pt_resp_all_w(3);

    vector<TH2D*> N_R_30(3);
    vector<TH2D*> N_R_60(3);
    vector<TH2D*> N_R_120(3);
    vector<TH2D*> N_R_240(3);
    vector<TH2D*> N_R_480(3);
    vector<TH2D*> N_R_960(3);
    vector<TH2D*> N_R_1920(3);
    vector<TH2D*> N_R_3840(3);

    vector<TH1D*> pt_resp_gaus(3);
    vector<TH1D*> pt_resp_nGenJetPF_w_gaus(3);
    vector<TH1D*> pt_resp_nGenJetPFSig_w_gaus(3);
    vector<TH2D*> bounds(3);

    vector<TH2D*> pt_resp_hist(3);
    vector<TH2D*> nGenPF_probs(3);
    vector<TH2D*> nGenPFSig_probs(3);
    vector<TH2D*> genJetMass(3);
    vector<TH2D*> jetGirth_probs(3);
    vector<TH2D*> genJetMass_probs(3);
    vector<TH2D*> genJetLHA_probs(3);
    vector<TH2D*> genJetWidth_probs(3);
    vector<TH2D*> genJetThrust_probs(3);
    vector<TH2D*> genJetMultiplicity_probs(3);

    vector<TH3D*> nGenJetPF_genJetWidth_probs(3);

    vector<TH1D*> dR_nPF_fromPV0(3);
    vector<TH1D*> dR_nPF_fromPV1(3);
    vector<TH1D*> dR_nPF_fromPV2(3);
    vector<TH1D*> dR_nPF_fromPV3(3);

    vector<TH2D*> neutralPtFrac(3);

    int type = 0;

   for (string str : {"gluon_", "quark_", "all_"}) {
    f->GetObject(char_add(str, "pt_resp"), pt_resp.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPF"), pt_resp_nGenJetPF.at(type));
    f->GetObject(char_add(str, "nGenPF_prof"), nGenPF.at(type));
    f->GetObject(char_add(str, "nGenPF_prof_w"), nGenPF_w.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPF_w"), pt_resp_nGenJetPF_w.at(type));
    f->GetObject(char_add(str, "genJetMass"), genJetMass.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPFSig_w"), pt_resp_nGenJetPFSig_w.at(type));
    f->GetObject(char_add(str, "pt_resp_jetGirth_w"), pt_resp_jetGirth_w.at(type));
    f2->GetObject(char_add(str, "nGenPFSig_prof"), nGenPFSig_prof.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetMass_w"), pt_resp_genJetMass_w.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetLHA_w"), pt_resp_genJetLHA_w.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetWidth_w"), pt_resp_genJetWidth_w.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetThrust_w"), pt_resp_genJetThrust_w.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetMultiplicity_w"), pt_resp_genJetMultiplicity_w.at(type));
    f->GetObject(char_add(str, "pt_resp_genJetPtD_w"), pt_resp_genJetPtD_w.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPF_genJetWidth_w"), pt_resp_nGenJetPF_genJetWidth_w.at(type));
    f->GetObject(char_add(str, "pt_resp_all_w"), pt_resp_all_w.at(type));

    f->GetObject(char_add(str, "N_R_30_norm"), N_R_30.at(type));
    f->GetObject(char_add(str, "N_R_60_norm"), N_R_60.at(type));
    f->GetObject(char_add(str, "N_R_120_norm"), N_R_120.at(type));
    f->GetObject(char_add(str, "N_R_240_norm"), N_R_240.at(type));
    f->GetObject(char_add(str, "N_R_480_norm"), N_R_480.at(type));
    f->GetObject(char_add(str, "N_R_960_norm"), N_R_960.at(type));
    f->GetObject(char_add(str, "N_R_1920_norm"), N_R_1920.at(type));
    f->GetObject(char_add(str, "N_R_3840_norm"), N_R_3840.at(type));

    f->GetObject(char_add(str, "pt_resp_gaus"), pt_resp_gaus.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPF_w_gaus"), pt_resp_nGenJetPF_w_gaus.at(type));
    f->GetObject(char_add(str, "pt_resp_nGenJetPFSig_w_gaus"), pt_resp_nGenJetPFSig_w_gaus.at(type));
    f->GetObject(char_add(str, "pt_resp_hist_norm"), pt_resp_hist.at(type));
    f->GetObject(char_add(str, "resp_gaus_bounds"), bounds.at(type));

    f2->GetObject(char_add(str, "nGenPF_probs"), nGenPF_probs.at(type));
    f2->GetObject(char_add(str, "nGenPFSig_probs"), nGenPFSig_probs.at(type));
    f2->GetObject(char_add(str, "jetGirth_probs"), jetGirth_probs.at(type));
    f2->GetObject(char_add(str, "genJetMass_probs"), genJetMass_probs.at(type));
    f2->GetObject(char_add(str, "genJet_LHA_probs"), genJetLHA_probs.at(type));
    f2->GetObject(char_add(str, "genJet_width_probs"), genJetWidth_probs.at(type));
    f2->GetObject(char_add(str, "genJet_thrust_probs"), genJetThrust_probs.at(type));
    f2->GetObject(char_add(str, "genJet_multiplicity_probs"), genJetMultiplicity_probs.at(type));
    f2->GetObject(char_add(str, "nGenJetPF_genJetWidth_probs"), nGenJetPF_genJetWidth_probs.at(type));

    f2->GetObject(char_add(str, "dR_nPF_fromPV0"), dR_nPF_fromPV0.at(type));
    f2->GetObject(char_add(str, "dR_nPF_fromPV1"), dR_nPF_fromPV1.at(type));
    f2->GetObject(char_add(str, "dR_nPF_fromPV2"), dR_nPF_fromPV2.at(type));
    f2->GetObject(char_add(str, "dR_nPF_fromPV3"), dR_nPF_fromPV3.at(type));

    f2->GetObject(char_add(str, "genJet_neutral_pT_fraction_norm"), neutralPtFrac.at(type));

    ++type;
   }

   vector<TH2D*> fnh_R_30(3);
   vector<TH2D*> fnh_R_60(3);
   vector<TH2D*> fnh_R_120(3);
   vector<TH2D*> fnh_R_240(3);
   vector<TH2D*> fnh_R_480(3);
   vector<TH2D*> fnh_R_960(3);
   vector<TH2D*> fnh_R_1920(3);
   vector<TH2D*> fnh_R_3840(3);

    type = G;
    for (string str : {"gluon_", "UD_", "S_"}) {
        f->GetObject(char_add(str, "fnh_R_30_norm"), fnh_R_30.at(type));
        f->GetObject(char_add(str, "fnh_R_60_norm"), fnh_R_60.at(type));
        f->GetObject(char_add(str, "fnh_R_120_norm"), fnh_R_120.at(type));
        f->GetObject(char_add(str, "fnh_R_240_norm"), fnh_R_240.at(type));
        f->GetObject(char_add(str, "fnh_R_480_norm"), fnh_R_480.at(type));
        f->GetObject(char_add(str, "fnh_R_960_norm"), fnh_R_960.at(type));
        f->GetObject(char_add(str, "fnh_R_1920_norm"), fnh_R_1920.at(type));
        f->GetObject(char_add(str, "fnh_R_3840_norm"), fnh_R_3840.at(type));
        ++type;
    }

   TH2D* jetNeutralPtFracUD;
   TH2D* jetNeutralPtFracS;
   f2->GetObject("UD_genJet_neutral_pT_fraction_norm", jetNeutralPtFracUD);
   f2->GetObject("S_genJet_neutral_pT_fraction_norm", jetNeutralPtFracS);

    TLatex *tex = new TLatex(); tex->SetNDC();
    tex->SetTextSize(0.05); tex->SetTextColor(kBlack);


    //Pt response
    TH1D* h = tdrHist("h", "Response", 0.9, 1.2,
    "gen p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c1 = tdrCanvas("c1", h, 4, 11, kSquare);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(pt_resp.at(ALL), "P", kOpenCircle, kBlack);
    tdrDraw(pt_resp_nGenJetPF_w.at(ALL), "", kOpenTriangleUp, kBlack);
//    tdrDraw(pt_resp_nGenJetPFSig_wg, "", kOpenTriangleDown, kBlack);
    tdrDraw(pt_resp.at(GLUON), "", kFullCircle, kRed);
    tdrDraw(pt_resp_nGenJetPF_w.at(GLUON), "", kFullTriangleUp, kRed);
//    tdrDraw(pt_resp.at(GLUON)_nGenJetPFSig_w, "", kFullTriangleDown, kOrange +7);
    tdrDraw(pt_resp.at(QUARK), "", kFullCircle, kBlue);
    tdrDraw(pt_resp_nGenJetPF_w.at(QUARK), "", kFullTriangleUp, kBlue);
//    tdrDraw(pt_resp_nGenJetPFSig_w.at(QUARK), "", kFullTriangleDown, kCyan +2);
    pt_resp.at(GLUON)->SetMarkerSize(1.2);
    pt_resp_nGenJetPF_w.at(GLUON)->SetMarkerSize(1.2);
    pt_resp.at(QUARK)->SetMarkerSize(1.2);
    pt_resp_nGenJetPF_w.at(QUARK)->SetMarkerSize(1.2);
    pt_resp_nGenJetPFSig_w.at(GLUON)->SetMarkerSize(1.2);
    pt_resp_nGenJetPFSig_w.at(QUARK)->SetMarkerSize(1.2);

    TLegend *leg1 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg1->AddEntry(pt_resp.at(ALL), "All", "PLE");
    leg1->AddEntry(pt_resp.at(GLUON), "Gluons", "PLE");
    leg1->AddEntry(pt_resp.at(QUARK), "Quarks", "PLE");
    leg1->AddEntry(pt_resp_nGenJetPF_w.at(GLUON), "Gluons weighted nGenJetPF", "PLE");
    leg1->AddEntry(pt_resp_nGenJetPF_w.at(QUARK), "Quarks weighted nGenJetPF", "PLE");
//    leg1->AddEntry(pt_resp_nGenJetPFSig_w.at(GLUON), "Gluons weighted (sig)", "PLE");
//    leg1->AddEntry(pt_resp_nGenJetPFSig_w.at(QUARK), "Quarks weighted (sig)", "PLE");

    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //Response gaus
    TCanvas* c7 = tdrCanvas("c7", h, 4, 11, kSquare);
    h->GetXaxis()->SetNoExponent();
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(pt_resp_gaus.at(ALL), "P", kOpenCircle, kBlack);
    tdrDraw(pt_resp_nGenJetPF_w_gaus.at(ALL), "P", kOpenTriangleUp, kBlack);
//    tdrDraw(pt_resp_nGenJetPFSig_wg_gaus, "P", kOpenTriangleDown, kBlack);
    tdrDraw(pt_resp_gaus.at(GLUON), "P", kFullCircle, kRed);
    tdrDraw(pt_resp_nGenJetPF_w_gaus.at(GLUON), "P", kFullTriangleUp, kRed);
//    tdrDraw(pt_resp_nGenJetPFSig_w_gaus.at(GLUON), "P", kFullTriangleDown, kOrange +7);
    tdrDraw(pt_resp_gaus.at(QUARK), "P", kFullCircle, kBlue);
    tdrDraw(pt_resp_nGenJetPF_w_gaus.at(QUARK), "P", kFullTriangleUp, kBlue);
//    tdrDraw(pt_resp_nGenJetPFSig_w_gaus.at(QUARK), "P", kFullTriangleDown, kCyan +2);
    pt_resp_gaus.at(GLUON)->SetMarkerSize(1.2);
    pt_resp_nGenJetPF_w_gaus.at(GLUON)->SetMarkerSize(1.2);
    pt_resp_gaus.at(QUARK)->SetMarkerSize(1.2);
    pt_resp_nGenJetPF_w_gaus.at(QUARK)->SetMarkerSize(1.2);
    pt_resp_nGenJetPFSig_w_gaus.at(GLUON)->SetMarkerSize(1.2);
    pt_resp_nGenJetPFSig_w_gaus.at(QUARK)->SetMarkerSize(1.2);

    TLegend *leg3 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg3->AddEntry(pt_resp_gaus.at(ALL), "All", "PLE");
    leg3->AddEntry(pt_resp_gaus.at(GLUON), "Gluons", "PLE");
    leg3->AddEntry(pt_resp_gaus.at(QUARK), "Quarks", "PLE");
    leg3->AddEntry(pt_resp_nGenJetPF_w_gaus.at(GLUON), "Gluons weighted", "PLE");
    leg3->AddEntry(pt_resp_nGenJetPF_w_gaus.at(QUARK), "Quarks weighted", "PLE");
//    leg3->AddEntry(pt_resp_nGenJetPFSig_w_gaus.at(GLUON), "Gluons weighted (sig)", "PLE");
//    leg3->AddEntry(pt_resp_nGenJetPFSig_w_gaus.at(QUARK), "Quarks weighted (sig)", "PLE");
 
    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //nGenJetPF profiles
    TH1D* h2 = tdrHist("h2", "nGenPF", 10, 70,
    "gen p_{T} (GeV)", 30, 3500);
    tdrCanvas("c2", h2, 4, 11,  kSquare);
    gPad->SetLogx();

    tdrDraw(nGenPF.at(GLUON),"", kFullCircle, kRed);
    tdrDraw(nGenPF.at(QUARK), "", kFullCircle, kBlue);
    tdrDraw(nGenPF_w.at(GLUON), "", kFullTriangleUp, kRed);
    tdrDraw(nGenPF_w.at(QUARK), "", kFullTriangleUp, kBlue);
    nGenPF.at(GLUON)->SetMarkerSize(1.5);
    nGenPF.at(QUARK)->SetMarkerSize(1.5);
    nGenPF_w.at(GLUON)->SetMarkerSize(1.5);
    nGenPF_w.at(QUARK)->SetMarkerSize(1.5);

    TLegend *leg2 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg2->AddEntry(nGenPF.at(GLUON), "Gluons", "PLE");
    leg2->AddEntry(nGenPF.at(QUARK), "Quarks", "PLE");
    leg2->AddEntry(nGenPF_w.at(GLUON), "Gluons, w", "PLE");
    leg2->AddEntry(nGenPF_w.at(QUARK), "Quarks, w", "PLE");

    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //Probablity distributions 3D
    TCanvas* c3 = new TCanvas("c3","c3",1400,700);
    c3->Divide(2,1);

    c3->cd(1);
    nGenPF_probs.at(GLUON)->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    nGenPF_probs.at(GLUON)->GetXaxis()->SetTitle("gen p{T} (GeV)");
    nGenPF_probs.at(GLUON)->GetYaxis()->SetTitle("nGenPF");
    nGenPF_probs.at(GLUON)->GetZaxis()->SetTitle("prob");
    nGenPF_probs.at(GLUON)->SetAxisRange(0,0.1,"Z");
    nGenPF_probs.at(GLUON)->SetTitle("gluons");
    nGenPF_probs.at(GLUON)->Draw("lego2");

    TText *gt = new TText(0.5, 0.5, "gluons");
    gt->Draw();

    c3->cd(2);
    nGenPF_probs.at(QUARK)->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    nGenPF_probs.at(QUARK)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenPF_probs.at(QUARK)->GetYaxis()->SetTitle("nGenPF");
    nGenPF_probs.at(QUARK)->GetZaxis()->SetTitle("P");
    nGenPF_probs.at(QUARK)->SetAxisRange(0,0.1, "Z");
    nGenPF_probs.at(QUARK)->SetTitle("quarks");
    nGenPF_probs.at(QUARK)->Draw("lego2");

    
    TText *qt = new TText(0.5, 0.5, "quarks");
    qt->Draw();

    //nGenJetPF probablity distributions 2D
    TCanvas* c4 = new TCanvas("c4","c4",1400,700);
    c4->Divide(2,1);

    c4->cd(1);
    nGenPF_probs.at(GLUON)->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    nGenPF_probs.at(GLUON)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenPF_probs.at(GLUON)->GetYaxis()->SetTitle("nGenPF");
    nGenPF_probs.at(GLUON)->GetZaxis()->SetTitle("prob");
    nGenPF_probs.at(GLUON)->SetAxisRange(0,0.1,"Z");
    nGenPF_probs.at(GLUON)->Draw("colz2");

    gt->Draw();

    c4->cd(2);
    nGenPF_probs.at(QUARK)->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    nGenPF_probs.at(QUARK)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenPF_probs.at(QUARK)->GetYaxis()->SetTitle("nGenPF");
    nGenPF_probs.at(QUARK)->GetZaxis()->SetTitle("prob");
    nGenPF_probs.at(QUARK)->SetAxisRange(0,0.1, "Z");
    nGenPF_probs.at(QUARK)->Draw("colz2");

    qt->Draw();

    //nGenJetPFSig probability distributions 2D
    TCanvas* c8 = new TCanvas("c8","c8",1400,700);
    c8->Divide(2,1);

    c8->cd(1);
    nGenPFSig_probs.at(GLUON)->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    nGenPFSig_probs.at(GLUON)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenPFSig_probs.at(GLUON)->GetYaxis()->SetTitle("nGenPF");
    nGenPFSig_probs.at(GLUON)->GetZaxis()->SetTitle("prob");
    nGenPFSig_probs.at(GLUON)->SetAxisRange(0,0.1,"Z");
    nGenPFSig_probs.at(GLUON)->Draw("colz");
    
    nGenPFSig_prof.at(GLUON)->SetAxisRange(30, 3500, "X");
    nGenPFSig_prof.at(GLUON)->SetMarkerStyle(kFullCircle);
    nGenPFSig_prof.at(GLUON)->Draw("same");

    TText *gt2 = new TText(0.5, 0.5, "gluons");
    gt2->Draw();

    c8->cd(2);
    nGenPFSig_probs.at(QUARK)->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    nGenPFSig_probs.at(QUARK)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenPFSig_probs.at(QUARK)->GetYaxis()->SetTitle("nGenPF");
    nGenPFSig_probs.at(QUARK)->GetZaxis()->SetTitle("prob");
    nGenPFSig_probs.at(QUARK)->SetAxisRange(0,0.1, "Z");
    nGenPFSig_probs.at(QUARK)->Draw("colz");

    nGenPFSig_prof.at(QUARK)->SetAxisRange(30, 3500, "X");
    nGenPFSig_prof.at(QUARK)->SetMarkerStyle(kFullCircle);
    nGenPFSig_prof.at(QUARK)->Draw("same");

    TText *qt2 = new TText(0.5, 0.5, "quarks");
    qt2->Draw();

    //nUEPF
    TCanvas* c9 = new TCanvas("c9", "c9", 700, 700);
    c9->SetRightMargin(0.15);
    nUEPF_nGenPF_hist->SetAxisRange(0, 100, "X");
    nUEPF_nGenPF_hist->Draw("colz");
    nUEPF_nGenPF_hist->GetXaxis()->SetTitle("nGenPF");
    nUEPF_nGenPF_hist->GetYaxis()->SetTitle("nUEPF_perA");

    TCanvas* c10 = new TCanvas("c10", "c10", 700, 700);
    c10->SetRightMargin(0.15);
    nUEPF_pt_hist->SetAxisRange(30, 3500, "X");
    nUEPF_pt_hist->Draw("colz");
    nUEPF_pt_hist->GetXaxis()->SetMoreLogLabels();
    nUEPF_pt_hist->GetXaxis()->SetNoExponent();
    nUEPF_pt_hist->GetXaxis()->SetTitle("gen p_{T}");
    nUEPF_pt_hist->GetYaxis()->SetTitle("nUEPF_perA");
    gPad->SetLogx();

    //Pt histograms
    TCanvas* c5 = new TCanvas("c5","c5",1400,700);
    c5->Divide(2,1);

    c5->cd(1);
    c5->cd(1)->SetRightMargin(0.15);
    pt_resp_hist.at(GLUON)->GetZaxis()->SetTitle("P");
    gPad->SetLogz();
    gPad->SetLogx();
    pt_resp_hist.at(GLUON)->SetAxisRange(30, 3500, "X");
    pt_resp_hist.at(GLUON)->Draw("colz");
    bounds.at(GLUON)->Draw("CONT SAME");
    pt_resp_hist.at(GLUON)->GetXaxis()->SetMoreLogLabels();
    pt_resp_hist.at(GLUON)->GetXaxis()->SetNoExponent();
//    pt_resp.at(GLUON)->SetMarkerColor(kBlack);
    pt_resp.at(GLUON)->Draw("SAME");
    pt_resp_gaus.at(GLUON)->Draw("P SAME");

    TLine* line_down =new TLine(30, 0.8, 3500, 0.8);
    line_down->SetLineColor(kYellow);
    line_down->SetLineWidth(2);

    c5->cd(2);
    c5->cd(2)->SetRightMargin(0.15);
    pt_resp_hist.at(QUARK)->GetZaxis()->SetTitle("P");
    gPad->SetLogz();
    gPad->SetLogx();
    pt_resp_hist.at(QUARK)->SetAxisRange(30, 3500, "X");
    pt_resp_hist.at(QUARK)->Draw("colz");
    bounds.at(QUARK)->Draw("CONT SAME");
    pt_resp_hist.at(QUARK)->GetXaxis()->SetMoreLogLabels();
    pt_resp_hist.at(QUARK)->GetXaxis()->SetNoExponent();
//    pt_resp.at(QUARK)->SetMarkerColor(kBlack);
    pt_resp.at(QUARK)->Draw("SAME");
    pt_resp_gaus.at(QUARK)->Draw("P SAME");

    // nPF_from_PV histograms
    Width_t w = 3;

    TH1D* h5 = tdrHist("h5", "N", 0, 25, "dR", 0, 1.5);
    tdrCanvas("c14", h5, 4, 11, true);

    dR_nPF_fromPV0.at(ALL)->SetLineColor(kBlack);
    dR_nPF_fromPV0.at(ALL)->SetLineWidth(w);
    dR_nPF_fromPV0.at(ALL)->Draw("same");
    dR_nPF_fromPV0.at(GLUON)->SetLineColor(kRed);
    dR_nPF_fromPV0.at(GLUON)->SetLineWidth(w);
    dR_nPF_fromPV0.at(GLUON)->Draw("same");
    dR_nPF_fromPV0.at(QUARK)->SetLineColor(kBlue);
    dR_nPF_fromPV0.at(QUARK)->SetLineWidth(w);
    dR_nPF_fromPV0.at(QUARK)->Draw("same");

    TLegend *leg4 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg4->AddEntry(dR_nPF_fromPV0.at(ALL), "All", "PLE");
    leg4->AddEntry(dR_nPF_fromPV0.at(GLUON), "Gluons", "PLE");
    leg4->AddEntry(dR_nPF_fromPV0.at(QUARK), "Quarks", "PLE");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 0");

    TH1D* h6 = tdrHist("h6", "N", 0, 10, "dR", 0, 1.5);
    tdrCanvas("c11", h6, 4, 11, true);

    dR_nPF_fromPV1.at(ALL)->SetLineColor(kBlack);
    dR_nPF_fromPV1.at(ALL)->SetLineWidth(w);
    dR_nPF_fromPV1.at(ALL)->Draw("same");
    dR_nPF_fromPV1.at(GLUON)->SetLineColor(kRed);
    dR_nPF_fromPV1.at(GLUON)->SetLineWidth(w);
    dR_nPF_fromPV1.at(GLUON)->Draw("same");
    dR_nPF_fromPV1.at(QUARK)->SetLineColor(kBlue);
    dR_nPF_fromPV1.at(QUARK)->SetLineWidth(w);
    dR_nPF_fromPV1.at(QUARK)->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 1");

    leg4->Draw();

    TH1D* h7 = tdrHist("h7", "N", 0, 3, "dR", 0, 1.5);
    tdrCanvas("c12", h7, 4, 11, true);
    
    dR_nPF_fromPV2.at(ALL)->SetLineColor(kBlack);
    dR_nPF_fromPV2.at(ALL)->SetLineWidth(w);
    dR_nPF_fromPV2.at(ALL)->Draw("same");
    dR_nPF_fromPV2.at(GLUON)->SetLineColor(kRed);
    dR_nPF_fromPV2.at(GLUON)->SetLineWidth(w);
    dR_nPF_fromPV2.at(GLUON)->Draw("same");
    dR_nPF_fromPV2.at(QUARK)->SetLineColor(kBlue);
    dR_nPF_fromPV2.at(QUARK)->SetLineWidth(w);
    dR_nPF_fromPV2.at(QUARK)->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 2");
    leg4->Draw();

    TH1D* h8 = tdrHist("h8", "N", 0, 25, "dR", 0, 1.5);
    tdrCanvas("c13", h8, 4, 11, true);

    dR_nPF_fromPV3.at(ALL)->SetLineColor(kBlack);
    dR_nPF_fromPV3.at(ALL)->SetLineWidth(w);
    dR_nPF_fromPV3.at(ALL)->Draw("same");
    dR_nPF_fromPV3.at(GLUON)->SetLineColor(kRed);
    dR_nPF_fromPV3.at(GLUON)->SetLineWidth(w);
    dR_nPF_fromPV3.at(GLUON)->Draw("same");
    dR_nPF_fromPV3.at(QUARK)->SetLineColor(kBlue);
    dR_nPF_fromPV3.at(QUARK)->SetLineWidth(w);
    dR_nPF_fromPV3.at(QUARK)->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 3");
    leg4->Draw();

    TH1D* h9 = tdrHist("h9", "N", 0, 10, "dR", 0, 1.5);
    tdrCanvas("c15", h9, 4, 11, true);
    gPad->SetLogy();

    dR_nPF_fromPV0.at(ALL)->SetLineColor(kBlack);
    dR_nPF_fromPV0.at(ALL)->Draw();
    dR_nPF_fromPV1.at(ALL)->SetLineColor(kMagenta);
    dR_nPF_fromPV1.at(ALL)->Draw("same");
    dR_nPF_fromPV2.at(ALL)->SetLineColor(kGreen);
    dR_nPF_fromPV2.at(ALL)->Draw("same");
    dR_nPF_fromPV3.at(ALL)->SetLineColor(kCyan+2);
    dR_nPF_fromPV3.at(ALL)->Draw("same");

    tex->DrawLatex(0.5, 0.95,"All");

/*
    TLegend *leg5 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg5->AddEntry(dR_nPF_fromPV0, "PF_{fromPV} = 0", "PLE");
    leg5->AddEntry(dR_nPF_fromPV1, "PF_{fromPV} = 1", "PLE");
    leg5->AddEntry(dR_nPF_fromPV2, "PF_{fromPV} = 2", "PLE");
    leg5->AddEntry(dR_nPF_fromPV3, "PF_{fromPV} = 3", "PLE");
    */

    vector<vector<TH2D*>> probs_histos = {genJetMass_probs, jetGirth_probs, genJetLHA_probs, genJetWidth_probs, genJetThrust_probs, genJetMultiplicity_probs};
    vector<const char*> y_label = {static_cast<const char*>("Mass"), static_cast<const char*>("Girth"), static_cast<const char*>("LHA"), static_cast<const char*>("Width"),
    static_cast<const char*>("Thrust"), static_cast<const char*>("Multiplicity")};
    for (int i = 0; i != probs_histos.size(); ++i) {
        TCanvas* c4 = new TCanvas(char_add("c", to_string(30 + i)), char_add("c", to_string(30 + i)),1400,700);
        c4->Divide(2,1);

        c4->cd(1);
        probs_histos.at(i).at(GLUON)->SetAxisRange(30, 3500, "X");
        gPad->SetLogx();
        probs_histos.at(i).at(GLUON)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
        probs_histos.at(i).at(GLUON)->GetYaxis()->SetTitle(y_label.at(i));
        probs_histos.at(i).at(GLUON)->GetZaxis()->SetTitle("prob");
        probs_histos.at(i).at(GLUON)->SetAxisRange(0,0.1,"Z");
        probs_histos.at(i).at(GLUON)->Draw("col");

        gt->Draw();

        c4->cd(2);
        probs_histos.at(i).at(QUARK)->SetAxisRange(30, 3500,"X");
        gPad->SetLogx();
        probs_histos.at(i).at(QUARK)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
        probs_histos.at(i).at(QUARK)->GetYaxis()->SetTitle(y_label.at(i));
        probs_histos.at(i).at(QUARK)->GetZaxis()->SetTitle("prob");
        probs_histos.at(i).at(QUARK)->SetAxisRange(0,0.1, "Z");
        probs_histos.at(i).at(QUARK)->Draw("col");
    }

    //Pt response weighted with different values
    vector<string> legend_str = {" mass", " girth", " LHA", " width", " thrust", " multiplicity", " pTD^2", " nGenJetPF+Width", " all"};
    vector<vector<TProfile*>> profiles = {pt_resp_genJetMass_w, pt_resp_jetGirth_w, pt_resp_genJetLHA_w, pt_resp_genJetWidth_w, pt_resp_genJetThrust_w, 
    pt_resp_genJetMultiplicity_w, pt_resp_genJetPtD_w, pt_resp_nGenJetPF_genJetWidth_w, pt_resp_all_w};
    vector<TProfile*> plotted_pt_prof;

    TH1D* allHist = pt_resp.at(ALL)->ProjectionX();
    pt_resp.at(GLUON)->Divide(allHist);
    pt_resp.at(QUARK)->Divide(allHist);
    pt_resp.at(ALL)->Divide(allHist);
    pt_resp.at(GLUON)->SetAxisRange(30, 3500, "X");
    pt_resp.at(QUARK)->SetAxisRange(30, 3500, "X");
    pt_resp.at(ALL)->SetAxisRange(30, 3500, "X");

    for (int i = 0; i != profiles.size(); ++i) {
        plotted_pt_prof = profiles.at(i);

        plotted_pt_prof.at(GLUON)->Divide(allHist);
        plotted_pt_prof.at(QUARK)->Divide(allHist);

        h->SetAxisRange(0.95, 1.1, "Y");
        TCanvas* c16 = tdrCanvas(char_add("c", to_string(16 + i)), h, 4, 11, kSquare);
        h->GetXaxis()->SetMoreLogLabels();
        gPad->SetLogx();

        plotted_pt_prof.at(GLUON)->SetAxisRange(30, 3500, "X");
        plotted_pt_prof.at(QUARK)->SetAxisRange(30, 3500, "X");

        tdrDraw(pt_resp.at(ALL), "P", kOpenCircle, kBlack);
        tdrDraw(pt_resp.at(GLUON), "P", kFullCircle, kRed);
        tdrDraw(plotted_pt_prof.at(GLUON), "P", kFullTriangleUp, kRed);
        tdrDraw(pt_resp.at(QUARK), "P", kFullCircle, kBlue);
        tdrDraw(plotted_pt_prof.at(QUARK), "P", kFullTriangleUp, kBlue);
        pt_resp.at(GLUON)->SetMarkerSize(1.2);
        plotted_pt_prof.at(GLUON)->SetMarkerSize(1.2);
        pt_resp.at(QUARK)->SetMarkerSize(1.2);
        plotted_pt_prof.at(QUARK)->SetMarkerSize(1.2);

        TLegend *leg5 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
        leg5->AddEntry(pt_resp.at(ALL), "All", "PLE");
        leg5->AddEntry(pt_resp.at(GLUON), "Gluons", "PLE");
        leg5->AddEntry(pt_resp.at(QUARK), "Quarks", "PLE");
        leg5->AddEntry(plotted_pt_prof.at(GLUON), char_add("Gluons w", legend_str.at(i)), "PLE");
        leg5->AddEntry(plotted_pt_prof.at(QUARK), char_add("Quarks w", legend_str.at(i)), "PLE");

        tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");
    }

    TCanvas* c40 = new TCanvas("c40", "c40", 1400, 700);
    c40->Divide(2,1);

    c40->cd(1);
    nGenJetPF_genJetWidth_probs.at(GLUON)->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    nGenJetPF_genJetWidth_probs.at(GLUON)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenJetPF_genJetWidth_probs.at(GLUON)->GetYaxis()->SetTitle("nGenJetPF");
    nGenJetPF_genJetWidth_probs.at(GLUON)->GetZaxis()->SetTitle("Width");
    //nGenJetPF_genJetWidth_probs.at(GLUON)->SetAxisRange(0,0.1,"Z");
    nGenJetPF_genJetWidth_probs.at(GLUON)->Draw("BOX2");

    gt->Draw();

    c40->cd(2);
    nGenJetPF_genJetWidth_probs.at(QUARK)->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    nGenJetPF_genJetWidth_probs.at(QUARK)->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenJetPF_genJetWidth_probs.at(QUARK)->GetYaxis()->SetTitle("nGenJetPF");
    nGenJetPF_genJetWidth_probs.at(QUARK)->GetZaxis()->SetTitle("Width");
    //nGenJetPF_genJetWidth_probs.at(QUARK)->SetAxisRange(0,0.1, "Z");
    nGenJetPF_genJetWidth_probs.at(QUARK)->Draw("BOX2");

    qt->Draw();
/*
    TFile* f3 = TFile::Open("outputGluonHistos_gen.root");

    vector<TH2D*> N_R_30_other(3);
    vector<TH2D*> N_R_60_other(3);
    vector<TH2D*> N_R_120_other(3);
    vector<TH2D*> N_R_240_other(3);
    vector<TH2D*> N_R_480_other(3);
    vector<TH2D*> N_R_960_other(3);
    vector<TH2D*> N_R_1920_other(3);
    vector<TH2D*> N_R_3840_other(3);

    type = 0;
    for (string str : {"gluon_", "quark_", "all_"}) {
        f3->GetObject(char_add(str, "N_R_30_norm"), N_R_30_other.at(type));
        f3->GetObject(char_add(str, "N_R_60_norm"), N_R_60_other.at(type));
        f3->GetObject(char_add(str, "N_R_120_norm"), N_R_120_other.at(type));
        f3->GetObject(char_add(str, "N_R_240_norm"), N_R_240_other.at(type));
        f3->GetObject(char_add(str, "N_R_480_norm"), N_R_480_other.at(type));
        f3->GetObject(char_add(str, "N_R_960_norm"), N_R_960_other.at(type));
        f3->GetObject(char_add(str, "N_R_1920_norm"), N_R_1920_other.at(type));
        f3->GetObject(char_add(str, "N_R_3840_norm"), N_R_3840_other.at(type));
        ++type;
    }

    Width_t width = 3;

    vector<TProfile*> profs_GQ(3);
    vector<TProfile*> profs_otherPt(3);
    vector<TH2D*> N_R_2D;
    vector<vector<TH2D*>> NR_histos = {N_R_30, N_R_60, N_R_120, N_R_240, N_R_480, N_R_960, N_R_1920, N_R_3840};
    vector<vector<TH2D*>> NR_histos_other = {N_R_30_other, N_R_60_other, N_R_120_other, N_R_240_other, N_R_480_other, N_R_960_other, N_R_1920_other, N_R_3840};
    for (int i = 0; i != NR_histos.size(); ++i) {
        N_R_2D = NR_histos.at(i);
        profs_GQ.at(GLUON) = N_R_2D.at(GLUON)->ProfileX(char_add("same_gluon", to_string(i)));
        profs_GQ.at(QUARK) = N_R_2D.at(QUARK)->ProfileX(char_add("same_quark", to_string(i)));
        profs_otherPt.at(GLUON) = NR_histos_other.at(i).at(GLUON)->ProfileX(char_add("other_gluon", to_string(i)));
        profs_otherPt.at(QUARK) = NR_histos_other.at(i).at(QUARK)->ProfileX(char_add("other_quark", to_string(i)));
        TCanvas* c50 = new TCanvas(char_add("c", to_string(50 + i)), char_add("c", to_string(50 + i)),1500,700);
        c50->Divide(2,1);

        c50->cd(1);
        N_R_2D.at(GLUON)->SetAxisRange(0.5, 1.8, "Y");
        N_R_2D.at(GLUON)->SetAxisRange(0, 130, "X");
        N_R_2D.at(GLUON)->GetXaxis()->SetTitle("N");
        N_R_2D.at(GLUON)->GetYaxis()->SetTitle("R");
        N_R_2D.at(GLUON)->GetZaxis()->SetTitle("count");
        N_R_2D.at(GLUON)->Draw("colz");
        profs_GQ.at(GLUON)->Draw("HISTPSAME");
        profs_GQ.at(QUARK)->Draw("HISTPSAME");
        profs_otherPt.at(GLUON)->Draw("HISTPSAME");

        gPad->SetLogz();

        tex->DrawLatex(0.15, 0.95, char_add("Gluons, correllation factor: ", to_string(N_R_2D.at(GLUON)->GetCorrelationFactor())));

        c50->cd(2);
        N_R_2D.at(QUARK)->SetAxisRange(0.5, 1.8, "Y");
        N_R_2D.at(QUARK)->SetAxisRange(0, 130,"X");
        N_R_2D.at(QUARK)->GetXaxis()->SetTitle("N");
        N_R_2D.at(QUARK)->GetYaxis()->SetTitle("R");
        N_R_2D.at(QUARK)->GetZaxis()->SetTitle("count");
        N_R_2D.at(QUARK)->Draw("colz");
        profs_GQ.at(GLUON)->Draw("HISTPSAME");
        profs_GQ.at(QUARK)->Draw("HISTPSAME");
        profs_otherPt.at(QUARK)->Draw("HISTP*SAME");

        gPad->SetLogz();

        profs_GQ.at(GLUON)->SetMarkerStyle(kFullCircle);
        profs_GQ.at(GLUON)->SetMarkerColor(kRed);
        profs_GQ.at(GLUON)->SetMarkerSize(0.8);
        profs_GQ.at(QUARK)->SetMarkerStyle(kFullCircle);
        profs_GQ.at(QUARK)->SetMarkerColor(kBlue);
        profs_GQ.at(QUARK)->SetMarkerSize(0.8);

        profs_otherPt.at(GLUON)->SetMarkerStyle(kFullCircle);
        profs_otherPt.at(GLUON)->SetMarkerColor(kBlack);
        profs_otherPt.at(GLUON)->SetMarkerSize(0.8);
        profs_otherPt.at(QUARK)->SetMarkerStyle(kFullCircle);
        profs_otherPt.at(QUARK)->SetMarkerColor(kBlack);
        profs_otherPt.at(QUARK)->SetMarkerSize(0.8);

        tex->DrawLatex(0.15, 0.95, char_add("Quarks, correllation factor: ", to_string(N_R_2D.at(QUARK)->GetCorrelationFactor())));
    }*/

    TCanvas* c60 = new TCanvas("c60", "60", 1500, 700);
    c60->Divide(2,1);

    c60->cd(1);
    neutralPtFrac.at(GLUON)->SetAxisRange(30, 3500, "X");
    neutralPtFrac.at(GLUON)->Draw("colz");
    neutralPtFrac.at(GLUON)->GetXaxis()->SetTitle("p_{T, jet}");
    neutralPtFrac.at(GLUON)->GetYaxis()->SetTitle("p_{T, neutral}/p_{T, jet}");
    gPad->SetLogx();
    gPad->SetLogz();

    c60->cd(2);
    neutralPtFrac.at(QUARK)->SetAxisRange(30, 3500, "X");
    neutralPtFrac.at(QUARK)->Draw("colz");
    neutralPtFrac.at(QUARK)->GetXaxis()->SetTitle("p_{T, jet}");
    neutralPtFrac.at(QUARK)->GetYaxis()->SetTitle("p_{T, neutral}/p_{T, jet}");
    gPad->SetLogx();
    gPad->SetLogz();

    TH1D* h11 = tdrHist("h11", "p_{T, neutral}/p_{T, jet}", 0, 1, "p_{T, jet} (GeV)" , 30, 3500);
    tdrCanvas("c61", h11, 4, 11, true);

    TProfile* neutralFrac_gluon_profile = neutralPtFrac.at(GLUON)->ProfileX("neutralFrac_gluon_profile");
    tdrDraw(neutralFrac_gluon_profile, "HISTP", kFullCircle, kRed);

    TProfile* neutralFrac_quark_profile = neutralPtFrac.at(QUARK)->ProfileX("neutralFrac_quark_profile");
    tdrDraw(neutralFrac_quark_profile, "HISTP", kFullCircle, kBlue);

    TProfile* neutralFrac_UD_profile = jetNeutralPtFracUD->ProfileX("neutralFrac_UD_profile");
    tdrDraw(neutralFrac_UD_profile, "HISTP", kFullCircle, kBlack);

    TProfile* neutralFrac_S_profile = jetNeutralPtFracS->ProfileX("neutralFrac_S_profile");
    tdrDraw(neutralFrac_S_profile, "HISTP", kFullCircle, kGreen);

    gPad->SetLogx();

    TLegend *leg5 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg5->AddEntry(neutralFrac_gluon_profile, "Gluons", "PLE");
    leg5->AddEntry(neutralFrac_quark_profile, "Quarks (UDS)", "PLE");
    leg5->AddEntry(neutralFrac_UD_profile, "UD", "PLE");
    leg5->AddEntry(neutralFrac_S_profile, "S", "PLE");


    vector<TProfile*> profs(3);
    vector<TH2D*> fnh_R_2D;
    vector<vector<TH2D*>> fnh_histos = {fnh_R_30, fnh_R_60, fnh_R_120, fnh_R_240, fnh_R_480, fnh_R_960, fnh_R_1920, fnh_R_3840};
    vector<string> pt_str = {"30 GeV", "60 GeV", "120 GeV", "240 GeV", "480 GeV", "960 GeV", "1920 GeV", "3840 GeV"};
    for (int i = 0; i != fnh_histos.size(); ++i) {
        fnh_R_2D = fnh_histos.at(i);
        profs.at(G) = fnh_R_2D.at(G)->ProfileX(char_add("gluon", to_string(i)));
        profs.at(UD) = fnh_R_2D.at(UD)->ProfileX(char_add("UD", to_string(i)));
        profs.at(S) = fnh_R_2D.at(S)->ProfileX(char_add("S", to_string(i)));
        TCanvas* c70 = new TCanvas(char_add("c", to_string(70 + i)), char_add("c", to_string(70 + i)),1600,900);
        c70->Divide(2,2);

        c70->cd(2);
        fnh_R_2D.at(GLUON)->SetAxisRange(0.5, 1.8, "Y");
        fnh_R_2D.at(GLUON)->SetAxisRange(0, 1, "X");
        fnh_R_2D.at(GLUON)->GetXaxis()->SetTitle("p_{T,neutral}/p_{T,jet}");
        fnh_R_2D.at(GLUON)->GetYaxis()->SetTitle("R");
        fnh_R_2D.at(GLUON)->GetZaxis()->SetTitle("count");
        fnh_R_2D.at(GLUON)->Draw("colz");
        profs.at(G)->Draw("HISTPSAME");
        profs.at(UD)->Draw("HISTPSAME");
        profs.at(S)->Draw("HISTPSAME");

        gPad->SetLogz();

        tex->DrawLatex(0.20, 0.95, "Gluons");

        c70->cd(3);
        fnh_R_2D.at(UD)->SetAxisRange(0.5, 1.8, "Y");
        fnh_R_2D.at(UD)->SetAxisRange(0, 1,"X");
        fnh_R_2D.at(UD)->GetXaxis()->SetTitle("p_{T,neutral}/p_{T,jet}");
        fnh_R_2D.at(UD)->GetYaxis()->SetTitle("R");
        fnh_R_2D.at(UD)->GetZaxis()->SetTitle("count");
        fnh_R_2D.at(UD)->Draw("colz");
        tex->DrawLatex(0.20, 0.95, "UD");
        profs.at(G)->Draw("HISTPSAME");
        profs.at(UD)->Draw("HISTPSAME");
        profs.at(S)->Draw("HISTPSAME");

        gPad->SetLogz();

        c70->cd(4);
        fnh_R_2D.at(S)->SetAxisRange(0.5, 1.8, "Y");
        fnh_R_2D.at(S)->SetAxisRange(0, 1,"X");
        fnh_R_2D.at(S)->GetXaxis()->SetTitle("p_{T,neutral}/p_{T,jet}");
        fnh_R_2D.at(S)->GetYaxis()->SetTitle("R");
        fnh_R_2D.at(S)->GetZaxis()->SetTitle("count");
        fnh_R_2D.at(S)->Draw("colz");
        tex->DrawLatex(0.20, 0.95, "S");
        profs.at(G)->Draw("HISTPSAME");
        profs.at(UD)->Draw("HISTPSAME");
        profs.at(S)->Draw("HISTPSAME");

        gPad->SetLogz();

        profs.at(GLUON)->SetMarkerStyle(kFullCircle);
        profs.at(GLUON)->SetMarkerColor(kRed);
        profs.at(GLUON)->SetMarkerSize(0.8);
        profs.at(QUARK)->SetMarkerStyle(kFullCircle);
        profs.at(QUARK)->SetMarkerColor(kBlue);
        profs.at(QUARK)->SetMarkerSize(0.8);
        profs.at(S)->SetMarkerStyle(kFullCircle);
        profs.at(S)->SetMarkerColor(kBlack);
        profs.at(S)->SetMarkerSize(0.8);

        c70->cd(1);
        TLatex *tex2 = new TLatex(); tex2->SetNDC();
        tex2->SetTextSize(0.1);
        tex2->DrawLatex(0.15, 0.8, pt_str.at(i).c_str());

        TLegend* leg70 = new TLegend(0.15,0.5,0.4,0.7);
        leg70->AddEntry(profs.at(GLUON), "Gluons", "PLE");
        leg70->AddEntry(profs.at(UD), "UD quarks", "PLE");
        leg70->AddEntry(profs.at(S), "S quarks", "PLE");
        leg70->Draw();
    }
}