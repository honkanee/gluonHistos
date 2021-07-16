#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TAxis.h>

#include "tdrstyle_mod15.C"

void drawGluonHistos() {
    setTDRStyle();

    TFile* f = TFile::Open("outputGluonHistos.root");
    TFile* f2 = TFile::Open("weightHistos.root");

    //All
    TProfile* pt_resp;
    TProfile* pt_resp_nGenJetPF_wg;
    TProfile* pt_resp_nGenJetPFSig_wg;

    TH1D* pt_resp_gaus;
    TH1D* pt_resp_nGenJetPF_wg_gaus;
    TH1D* pt_resp_nGenJetPFSig_wg_gaus;
    TH2D* nUEPF_pt_hist;
    TH2D* nUEPF_nGenPF_hist;

    TH1D* dR_nPF_fromPV0;
    TH1D* dR_nPF_fromPV1;
    TH1D* dR_nPF_fromPV2;
    TH1D* dR_nPF_fromPV3;

    f->GetObject("pt_resp", pt_resp);
    f->GetObject("pt_resp_nGenJetPF_wg", pt_resp_nGenJetPF_wg);
    f->GetObject("pt_resp_nGenJetPFSig_wg", pt_resp_nGenJetPFSig_wg);

    f->GetObject("pt_resp_gaus", pt_resp_gaus);
    f->GetObject("pt_resp_nGenJetPF_wg_gaus", pt_resp_nGenJetPF_wg_gaus);
    f->GetObject("pt_resp_nGenJetPFSig_wg_gaus", pt_resp_nGenJetPFSig_wg_gaus);

    f2->GetObject("nUEPF_perA_genJetPt_probs", nUEPF_pt_hist);
    f2->GetObject("nUEPF_perA_nGenPF_probs", nUEPF_nGenPF_hist);

    f2->GetObject("dR_nPF_fromPV0", dR_nPF_fromPV0);
    f2->GetObject("dR_nPF_fromPV1", dR_nPF_fromPV1);
    f2->GetObject("dR_nPF_fromPV2", dR_nPF_fromPV2);
    f2->GetObject("dR_nPF_fromPV3", dR_nPF_fromPV3);

    // Gluons
    TProfile* gluon_pt_resp;
    TProfile* gluon_pt_resp_nGenJetPF;
    TProfile* gluon_nGenPF;
    TProfile* gluon_nGenPF_w;
    TProfile* gluon_pt_resp_nGenJetPF_w;
    TProfile* gluon_pt_resp_nGenJetPFSig_w;
    TProfile* gluon_nGenPFSig_prof;
    TProfile* gluon_pt_resp_jetGirth_w;

    TH1D* gluon_pt_resp_gaus;
    TH1D* gluon_pt_resp_nGenJetPF_w_gaus;
    TH1D* gluon_pt_resp_nGenJetPFSig_w_gaus;
    TH2D* gbounds;

    TH2D* gluon_pt_resp_hist;
    TH2D* gluon_nGenPF_probs;
    TH2D* gluon_nGenPFSig_probs;
    TH2D* gluon_genJetMass;

    TH1D* gluon_dR_nPF_fromPV0;
    TH1D* gluon_dR_nPF_fromPV1;
    TH1D* gluon_dR_nPF_fromPV2;
    TH1D* gluon_dR_nPF_fromPV3;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF_prof", gluon_nGenPF);
    f->GetObject("gluon_nGenPF_prof_w", gluon_nGenPF_w);
    f->GetObject("gluon_pt_resp_nGenJetPF_w", gluon_pt_resp_nGenJetPF_w);
    f->GetObject("gluon_genJetMass", gluon_genJetMass);
    f->GetObject("gluon_pt_resp_nGenJetPFSig_w", gluon_pt_resp_nGenJetPFSig_w);
    f->GetObject("gluon_pt_resp_jetGirth_w", gluon_pt_resp_jetGirth_w);
    f2->GetObject("gluon_nGenPFSig_prof", gluon_nGenPFSig_prof);

    f->GetObject("gluon_pt_resp_gaus", gluon_pt_resp_gaus);
    f->GetObject("gluon_pt_resp_nGenJetPF_w_gaus", gluon_pt_resp_nGenJetPF_w_gaus);
    f->GetObject("gluon_pt_resp_nGenJetPFSig_w_gaus", gluon_pt_resp_nGenJetPFSig_w_gaus);
    f->GetObject("gluon_pt_resp_hist_norm", gluon_pt_resp_hist);
    f->GetObject("gluon_resp_gaus_bounds", gbounds);

    f2->GetObject("gluon_nGenPF_probs", gluon_nGenPF_probs);
    f2->GetObject("gluon_nGenPFSig_probs", gluon_nGenPFSig_probs);

    f2->GetObject("gluon_dR_nPF_fromPV0", gluon_dR_nPF_fromPV0);
    f2->GetObject("gluon_dR_nPF_fromPV1", gluon_dR_nPF_fromPV1);
    f2->GetObject("gluon_dR_nPF_fromPV2", gluon_dR_nPF_fromPV2);
    f2->GetObject("gluon_dR_nPF_fromPV3", gluon_dR_nPF_fromPV3);

    // Quarks
    TProfile* quark_pt_resp;
    TProfile* quark_pt_resp_nGenJetPF;
    TProfile* quark_nGenPF;
    TProfile* quark_nGenPF_w;
    TProfile* quark_pt_resp_nGenJetPF_w;
    TProfile* quark_pt_resp_nGenJetPFSig_w;
    TProfile* quark_nGenPFSig_prof;
    TProfile* quark_pt_resp_jetGirth_w;

    TH1D* quark_pt_resp_gaus;
    TH1D* quark_pt_resp_nGenJetPF_w_gaus;
    TH1D* quark_pt_resp_nGenJetPFSig_w_gaus;
    TH2D* qbounds;

    TH2D* quark_pt_resp_hist;
    TH2D* quark_nGenPF_probs;
    TH2D* quark_nGenPFSig_probs;
    TH2D* quark_genJetMass;
    TH2D* nGenJetPF_likelyhood;

    TH1D* quark_dR_nPF_fromPV0;
    TH1D* quark_dR_nPF_fromPV1;
    TH1D* quark_dR_nPF_fromPV2;
    TH1D* quark_dR_nPF_fromPV3;

    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF_prof", quark_nGenPF);
    f->GetObject("quark_nGenPF_prof_w", quark_nGenPF_w);
    f->GetObject("quark_pt_resp_nGenJetPF_w", quark_pt_resp_nGenJetPF_w);
    f->GetObject("quark_genJetMass", quark_genJetMass);
    f->GetObject("quark_pt_resp_nGenJetPFSig_w", quark_pt_resp_nGenJetPFSig_w);
    f->GetObject("quark_pt_resp_jetGirth_w", quark_pt_resp_jetGirth_w);
    f2->GetObject("quark_nGenPFSig_prof", quark_nGenPFSig_prof);

    f->GetObject("quark_pt_resp_gaus", quark_pt_resp_gaus);
    f->GetObject("quark_pt_resp_nGenJetPF_w_gaus", quark_pt_resp_nGenJetPF_w_gaus);
    f->GetObject("quark_pt_resp_nGenJetPFSig_w_gaus", quark_pt_resp_nGenJetPFSig_w_gaus);
    f->GetObject("quark_pt_resp_hist_norm", quark_pt_resp_hist);
    f->GetObject("quark_resp_gaus_bounds", qbounds);

    f2->GetObject("quark_nGenPF_probs", quark_nGenPF_probs);
    f2->GetObject("quark_nGenPFSig_probs", quark_nGenPFSig_probs);

    f2->GetObject("quark_dR_nPF_fromPV0", quark_dR_nPF_fromPV0);
    f2->GetObject("quark_dR_nPF_fromPV1", quark_dR_nPF_fromPV1);
    f2->GetObject("quark_dR_nPF_fromPV2", quark_dR_nPF_fromPV2);
    f2->GetObject("quark_dR_nPF_fromPV3", quark_dR_nPF_fromPV3);

    f->GetObject("nGenJetPF_likelyhood", nGenJetPF_likelyhood);

    //Pt response
    TH1D* h = tdrHist("h", "Response", 0.9, 1.2,
    "gen p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c1 = tdrCanvas("c1", h, 4, 11, kSquare);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(pt_resp, "P", kOpenCircle, kBlack);
    tdrDraw(pt_resp_nGenJetPF_wg, "", kOpenTriangleUp, kBlack);
//    tdrDraw(pt_resp_nGenJetPFSig_wg, "", kOpenTriangleDown, kBlack);
    tdrDraw(gluon_pt_resp, "", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kRed);
//    tdrDraw(gluon_pt_resp_nGenJetPFSig_w, "", kFullTriangleDown, kOrange +7);
    tdrDraw(quark_pt_resp, "", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kBlue);
//    tdrDraw(quark_pt_resp_nGenJetPFSig_w, "", kFullTriangleDown, kCyan +2);
    gluon_pt_resp->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);
    quark_pt_resp->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPFSig_w->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPFSig_w->SetMarkerSize(1.2);

    TLegend *leg1 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg1->AddEntry(pt_resp, "All", "PLE");
    leg1->AddEntry(gluon_pt_resp, "Gluons", "PLE");
    leg1->AddEntry(quark_pt_resp, "Quarks", "PLE");
    leg1->AddEntry(gluon_pt_resp_nGenJetPF_w, "Gluons weighted", "PLE");
    leg1->AddEntry(quark_pt_resp_nGenJetPF_w, "Quarks weighted", "PLE");
//    leg1->AddEntry(gluon_pt_resp_nGenJetPFSig_w, "Gluons weighted (sig)", "PLE");
//    leg1->AddEntry(quark_pt_resp_nGenJetPFSig_w, "Quarks weighted (sig)", "PLE");

    TLatex *tex = new TLatex(); tex->SetNDC();
    tex->SetTextSize(0.05); tex->SetTextColor(kBlack);
    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //Response gaus
    TCanvas* c7 = tdrCanvas("c7", h, 4, 11, kSquare);
    h->GetXaxis()->SetNoExponent();
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(pt_resp_gaus, "P", kOpenCircle, kBlack);
    tdrDraw(pt_resp_nGenJetPF_wg_gaus, "P", kOpenTriangleUp, kBlack);
//    tdrDraw(pt_resp_nGenJetPFSig_wg_gaus, "P", kOpenTriangleDown, kBlack);
    tdrDraw(gluon_pt_resp_gaus, "P", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPF_w_gaus, "P", kFullTriangleUp, kRed);
//    tdrDraw(gluon_pt_resp_nGenJetPFSig_w_gaus, "P", kFullTriangleDown, kOrange +7);
    tdrDraw(quark_pt_resp_gaus, "P", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPF_w_gaus, "P", kFullTriangleUp, kBlue);
//    tdrDraw(quark_pt_resp_nGenJetPFSig_w_gaus, "P", kFullTriangleDown, kCyan +2);
    gluon_pt_resp_gaus->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPF_w_gaus->SetMarkerSize(1.2);
    quark_pt_resp_gaus->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPF_w_gaus->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPFSig_w_gaus->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPFSig_w_gaus->SetMarkerSize(1.2);

    TLegend *leg3 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg3->AddEntry(pt_resp_gaus, "All", "PLE");
    leg3->AddEntry(gluon_pt_resp_gaus, "Gluons", "PLE");
    leg3->AddEntry(quark_pt_resp_gaus, "Quarks", "PLE");
    leg3->AddEntry(gluon_pt_resp_nGenJetPF_w_gaus, "Gluons weighted", "PLE");
    leg3->AddEntry(quark_pt_resp_nGenJetPF_w_gaus, "Quarks weighted", "PLE");
//    leg3->AddEntry(gluon_pt_resp_nGenJetPFSig_w_gaus, "Gluons weighted (sig)", "PLE");
//    leg3->AddEntry(quark_pt_resp_nGenJetPFSig_w_gaus, "Quarks weighted (sig)", "PLE");
 
    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //nGenJetPF profiles
    TH1D* h2 = tdrHist("h2", "nGenPF", 10, 70,
    "gen p_{T} (GeV)", 30, 3500);
    tdrCanvas("c2", h2, 4, 11,  kSquare);
    gPad->SetLogx();

    tdrDraw(gluon_nGenPF,"", kFullCircle, kRed);
    tdrDraw(quark_nGenPF, "", kFullCircle, kBlue);
    tdrDraw(gluon_nGenPF_w, "", kFullTriangleUp, kRed);
    tdrDraw(quark_nGenPF_w, "", kFullTriangleUp, kBlue);
    gluon_nGenPF->SetMarkerSize(1.5);
    quark_nGenPF->SetMarkerSize(1.5);
    gluon_nGenPF_w->SetMarkerSize(1.5);
    quark_nGenPF_w->SetMarkerSize(1.5);

    TLegend *leg2 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg2->AddEntry(gluon_nGenPF, "Gluons", "PLE");
    leg2->AddEntry(quark_nGenPF, "Quarks", "PLE");
    leg2->AddEntry(gluon_nGenPF_w, "Gluons, w", "PLE");
    leg2->AddEntry(quark_nGenPF_w, "Quarks, w", "PLE");

    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //Probablity distributions 3D
    TCanvas* c3 = new TCanvas("c3","c3",1400,700);
    c3->Divide(2,1);

    c3->cd(1);
    gluon_nGenPF_probs->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    gluon_nGenPF_probs->GetXaxis()->SetTitle("gen p{T} (GeV)");
    gluon_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    gluon_nGenPF_probs->GetZaxis()->SetTitle("prob");
    gluon_nGenPF_probs->SetAxisRange(0,0.1,"Z");
    gluon_nGenPF_probs->SetTitle("gluons");
    gluon_nGenPF_probs->Draw("lego2");

    TText *gt = new TText(0.5, 0.5, "gluons");
    gt->Draw();

    c3->cd(2);
    quark_nGenPF_probs->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    quark_nGenPF_probs->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    quark_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPF_probs->GetZaxis()->SetTitle("P");
    quark_nGenPF_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPF_probs->SetTitle("quarks");
    quark_nGenPF_probs->Draw("lego2");

    
    TText *qt = new TText(0.5, 0.5, "quarks");
    qt->Draw();

    //nGenJetPF probablity distributions 2D
    TCanvas* c4 = new TCanvas("c4","c4",1400,700);
    c4->Divide(2,1);

    c4->cd(1);
    gluon_nGenPF_probs->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    gluon_nGenPF_probs->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    gluon_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    gluon_nGenPF_probs->GetZaxis()->SetTitle("prob");
    gluon_nGenPF_probs->SetAxisRange(0,0.1,"Z");
    gluon_nGenPF_probs->Draw("colz2");

    gt->Draw();

    c4->cd(2);
    quark_nGenPF_probs->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    quark_nGenPF_probs->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    quark_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPF_probs->GetZaxis()->SetTitle("prob");
    quark_nGenPF_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPF_probs->Draw("colz2");

    qt->Draw();

    //nGenJetPFSig probability distributions 2D
    TCanvas* c8 = new TCanvas("c8","c8",1400,700);
    c8->Divide(2,1);

    c8->cd(1);
    gluon_nGenPFSig_probs->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    gluon_nGenPFSig_probs->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    gluon_nGenPFSig_probs->GetYaxis()->SetTitle("nGenPF");
    gluon_nGenPFSig_probs->GetZaxis()->SetTitle("prob");
    gluon_nGenPFSig_probs->SetAxisRange(0,0.1,"Z");
    gluon_nGenPFSig_probs->Draw("colz");
    
    gluon_nGenPFSig_prof->SetAxisRange(30, 3500, "X");
    gluon_nGenPFSig_prof->SetMarkerStyle(kFullCircle);
    gluon_nGenPFSig_prof->Draw("same");

    TText *gt2 = new TText(0.5, 0.5, "gluons");
    gt2->Draw();

    c8->cd(2);
    quark_nGenPFSig_probs->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    quark_nGenPFSig_probs->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    quark_nGenPFSig_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPFSig_probs->GetZaxis()->SetTitle("prob");
    quark_nGenPFSig_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPFSig_probs->Draw("colz");

    quark_nGenPFSig_prof->SetAxisRange(30, 3500, "X");
    quark_nGenPFSig_prof->SetMarkerStyle(kFullCircle);
    quark_nGenPFSig_prof->Draw("same");

    TText *qt2 = new TText(0.5, 0.5, "quarks");
    qt2->Draw();

    //likelihood
    TH1D* h3 = tdrHist("h3", "nGenPF", 10, 70,
    "p_{T} (GeV)", 30, 3500);
    h3->GetXaxis()->SetNoExponent();
    h3->GetXaxis()->SetMoreLogLabels();
    tdrCanvas("c6", h3, 4, 11,  kSquare);
    gPad->SetLogx();

    nGenJetPF_likelyhood->SetAxisRange(30, 3500, "X");
    nGenJetPF_likelyhood->SetAxisRange(0, 150, "Y");
    nGenJetPF_likelyhood->GetXaxis()->SetTitle("gen p_{T} (GeV)");
    nGenJetPF_likelyhood->GetYaxis()->SetTitle("n_{GenJetPF}");
    nGenJetPF_likelyhood->Draw("colz");

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
    gluon_pt_resp_hist->GetZaxis()->SetTitle("P");
    gPad->SetLogz();
    gPad->SetLogx();
    gluon_pt_resp_hist->SetAxisRange(30, 3500, "X");
    gluon_pt_resp_hist->Draw("colz");
    gbounds->Draw("CONT SAME");
    gluon_pt_resp_hist->GetXaxis()->SetMoreLogLabels();
    gluon_pt_resp_hist->GetXaxis()->SetNoExponent();
//    gluon_pt_resp->SetMarkerColor(kBlack);
    gluon_pt_resp->Draw("SAME");
    gluon_pt_resp_gaus->Draw("P SAME");

    TLine* line_down =new TLine(30, 0.8, 3500, 0.8);
    line_down->SetLineColor(kYellow);
    line_down->SetLineWidth(2);

    c5->cd(2);
    c5->cd(2)->SetRightMargin(0.15);
    quark_pt_resp_hist->GetZaxis()->SetTitle("P");
    gPad->SetLogz();
    gPad->SetLogx();
    quark_pt_resp_hist->SetAxisRange(30, 3500, "X");
    quark_pt_resp_hist->Draw("colz");
    qbounds->Draw("CONT SAME");
    quark_pt_resp_hist->GetXaxis()->SetMoreLogLabels();
    quark_pt_resp_hist->GetXaxis()->SetNoExponent();
//    quark_pt_resp->SetMarkerColor(kBlack);
    quark_pt_resp->Draw("SAME");
    quark_pt_resp_gaus->Draw("P SAME");

    // nPF_from_PV histograms
    Width_t w = 3;

    TH1D* h5 = tdrHist("h5", "N", 0, 25, "dR", 0, 1.5);
    tdrCanvas("c14", h5, 4, 11, true);

    dR_nPF_fromPV0->SetLineColor(kBlack);
    dR_nPF_fromPV0->SetLineWidth(w);
    dR_nPF_fromPV0->Draw("same");
    gluon_dR_nPF_fromPV0->SetLineColor(kRed);
    gluon_dR_nPF_fromPV0->SetLineWidth(w);
    gluon_dR_nPF_fromPV0->Draw("same");
    quark_dR_nPF_fromPV0->SetLineColor(kBlue);
    quark_dR_nPF_fromPV0->SetLineWidth(w);
    quark_dR_nPF_fromPV0->Draw("same");

    TLegend *leg4 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg4->AddEntry(dR_nPF_fromPV0, "All", "PLE");
    leg4->AddEntry(gluon_dR_nPF_fromPV0, "Gluons", "PLE");
    leg4->AddEntry(quark_dR_nPF_fromPV0, "Quarks", "PLE");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 0");

    TH1D* h6 = tdrHist("h6", "N", 0, 10, "dR", 0, 1.5);
    tdrCanvas("c11", h6, 4, 11, true);

    dR_nPF_fromPV1->SetLineColor(kBlack);
    dR_nPF_fromPV1->SetLineWidth(w);
    dR_nPF_fromPV1->Draw("same");
    gluon_dR_nPF_fromPV1->SetLineColor(kRed);
    gluon_dR_nPF_fromPV1->SetLineWidth(w);
    gluon_dR_nPF_fromPV1->Draw("same");
    quark_dR_nPF_fromPV1->SetLineColor(kBlue);
    quark_dR_nPF_fromPV1->SetLineWidth(w);
    quark_dR_nPF_fromPV1->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 1");

    leg4->Draw();

    TH1D* h7 = tdrHist("h7", "N", 0, 3, "dR", 0, 1.5);
    tdrCanvas("c12", h7, 4, 11, true);
    
    dR_nPF_fromPV2->SetLineColor(kBlack);
    dR_nPF_fromPV2->SetLineWidth(w);
    dR_nPF_fromPV2->Draw("same");
    gluon_dR_nPF_fromPV2->SetLineColor(kRed);
    gluon_dR_nPF_fromPV2->SetLineWidth(w);
    gluon_dR_nPF_fromPV2->Draw("same");
    quark_dR_nPF_fromPV2->SetLineColor(kBlue);
    quark_dR_nPF_fromPV2->SetLineWidth(w);
    quark_dR_nPF_fromPV2->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 2");
    leg4->Draw();

    TH1D* h8 = tdrHist("h8", "N", 0, 25, "dR", 0, 1.5);
    tdrCanvas("c13", h8, 4, 11, true);

    dR_nPF_fromPV3->SetLineColor(kBlack);
    dR_nPF_fromPV3->SetLineWidth(w);
    dR_nPF_fromPV3->Draw("same");
    gluon_dR_nPF_fromPV3->SetLineColor(kRed);
    gluon_dR_nPF_fromPV3->SetLineWidth(w);
    gluon_dR_nPF_fromPV3->Draw("same");
    quark_dR_nPF_fromPV3->SetLineColor(kBlue);
    quark_dR_nPF_fromPV3->SetLineWidth(w);
    quark_dR_nPF_fromPV3->Draw("same");

    tex->DrawLatex(0.5, 0.95,"PF_{fromPV} = 3");
    leg4->Draw();

/*
    TH1D* h9 = tdrHist("h9", "N", 0, 10, "dR", 0, 1.5);
    tdrCanvas("c15", h9, 4, 11, true);
    gPad->SetLogy();

    dR_nPF_fromPV0->SetLineColor(kBlack);
    dR_nPF_fromPV0->Draw();
    dR_nPF_fromPV1->SetLineColor(kMagenta);
    dR_nPF_fromPV1->Draw("same");
    dR_nPF_fromPV2->SetLineColor(kGreen);
    dR_nPF_fromPV2->Draw("same");
    dR_nPF_fromPV3->SetLineColor(kCyan+2);
    dR_nPF_fromPV3->Draw("same");

    tex->DrawLatex(0.5, 0.95,"All");

    TLegend *leg5 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg5->AddEntry(dR_nPF_fromPV0, "PF_{fromPV} = 0", "PLE");
    leg5->AddEntry(dR_nPF_fromPV1, "PF_{fromPV} = 1", "PLE");
    leg5->AddEntry(dR_nPF_fromPV2, "PF_{fromPV} = 2", "PLE");
    leg5->AddEntry(dR_nPF_fromPV3, "PF_{fromPV} = 3", "PLE");
    */

    //Pt response jetGirth weighted
    TCanvas* c16 = tdrCanvas("c16", h, 4, 11, kSquare);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(pt_resp, "P", kOpenCircle, kBlack);
    tdrDraw(gluon_pt_resp, "", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_jetGirth_w, "", kFullTriangleUp, kRed);
    tdrDraw(quark_pt_resp, "", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_jetGirth_w, "", kFullTriangleUp, kBlue);
    gluon_pt_resp->SetMarkerSize(1.2);
    gluon_pt_resp_jetGirth_w->SetMarkerSize(1.2);
    quark_pt_resp->SetMarkerSize(1.2);
    quark_pt_resp_jetGirth_w->SetMarkerSize(1.2);

    TLegend *leg5 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg5->AddEntry(pt_resp, "All", "PLE");
    leg5->AddEntry(gluon_pt_resp, "Gluons", "PLE");
    leg5->AddEntry(quark_pt_resp, "Quarks", "PLE");
    leg5->AddEntry(gluon_pt_resp_jetGirth_w, "Gluons weighted (jetGirth)", "PLE");
    leg5->AddEntry(quark_pt_resp_jetGirth_w, "Quarks weighted (jetGirth)", "PLE");

    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");
}