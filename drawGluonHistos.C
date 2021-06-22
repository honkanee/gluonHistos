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

    // Gluons
    TProfile* gluon_pt_resp;
    TProfile* gluon_pt_resp_nGenJetPF;
    TProfile* gluon_nGenPF;
    TProfile* gluon_nGenPF_w;
    TProfile* gluon_pt_resp_nGenJetPF_w;
    TProfile* gluon_pt_resp_nGenJetPFSig_w;
    TProfile* gluon_nGenPFSig_prof;

    TH1D* gluon_pt_resp_gaus;
    TH1D* gluon_pt_resp_nGenJetPF_w_gaus;
    TH1D* gluon_pt_resp_nGenJetPFSig_w_gaus;

    TH2D* gluon_nGenPF_probs;
    TH2D* gluon_nGenPFSig_probs;
    TH2D* gluon_genJetMass;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF_prof", gluon_nGenPF);
    f->GetObject("gluon_nGenPF_prof_w", gluon_nGenPF_w);
    f->GetObject("gluon_pt_resp_nGenJetPF_w", gluon_pt_resp_nGenJetPF_w);
    f->GetObject("gluon_genJetMass", gluon_genJetMass);
    f->GetObject("gluon_pt_resp_nGenJetPFSig_w", gluon_pt_resp_nGenJetPFSig_w);
    f2->GetObject("gluon_nGenPFSig_prof", gluon_nGenPFSig_prof);

    f->GetObject("gluon_pt_resp_hist_1", gluon_pt_resp_gaus);
    f->GetObject("gluon_pt_resp_nGenJetPF_w_hist_1", gluon_pt_resp_nGenJetPF_w_gaus);
    f->GetObject("gluon_pt_resp_nGenJetPFSig_w_hist_1", gluon_pt_resp_nGenJetPFSig_w_gaus);

    f2->GetObject("gluon_nGenPF_probs", gluon_nGenPF_probs);
    f2->GetObject("gluon_nGenPFSig_probs", gluon_nGenPFSig_probs);


    // Quarks
    TProfile* quark_pt_resp;
    TProfile* quark_pt_resp_nGenJetPF;
    TProfile* quark_nGenPF;
    TProfile* quark_nGenPF_w;
    TProfile* quark_pt_resp_nGenJetPF_w;
    TProfile* quark_pt_resp_nGenJetPFSig_w;
    TProfile* quark_nGenPFSig_prof;

    TH1D* quark_pt_resp_gaus;
    TH1D* quark_pt_resp_nGenJetPF_w_gaus;
    TH1D* quark_pt_resp_nGenJetPFSig_w_gaus;

    TH2D* quark_nGenPF_probs;
    TH2D* quark_nGenPFSig_probs;
    TH2D* quark_genJetMass;
    TH2D* nGenJetPF_likelyhood;

    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF_prof", quark_nGenPF);
    f->GetObject("quark_nGenPF_prof_w", quark_nGenPF_w);
    f->GetObject("quark_pt_resp_nGenJetPF_w", quark_pt_resp_nGenJetPF_w);
    f->GetObject("quark_genJetMass", quark_genJetMass);
    f->GetObject("quark_pt_resp_nGenJetPFSig_w", quark_pt_resp_nGenJetPFSig_w);
    f2->GetObject("quark_nGenPFSig_prof", quark_nGenPFSig_prof);

    f->GetObject("quark_pt_resp_hist_1", quark_pt_resp_gaus);
    f->GetObject("quark_pt_resp_nGenJetPF_w_hist_1", quark_pt_resp_nGenJetPF_w_gaus);
    f->GetObject("quark_pt_resp_nGenJetPFSig_w_hist_1", quark_pt_resp_nGenJetPFSig_w_gaus);

    f2->GetObject("quark_nGenPF_probs", quark_nGenPF_probs);
    f2->GetObject("quark_nGenPFSig_probs", quark_nGenPFSig_probs);

    f->GetObject("nGenJetPF_likelyhood", nGenJetPF_likelyhood);

    //Pt response
    TH1D* h = tdrHist("h", "Response", 0.9, 1.2,
    "gen p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c1 = tdrCanvas("c1", h, 4, 11, kSquare);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(gluon_pt_resp, "P", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPFSig_w, "", kFullTriangleDown, kOrange +7);
    tdrDraw(quark_pt_resp, "", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPFSig_w, "", kFullTriangleDown, kCyan +2);
    gluon_pt_resp->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);
    quark_pt_resp->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPFSig_w->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPFSig_w->SetMarkerSize(1.2);

    TLegend *leg1 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg1->AddEntry(gluon_pt_resp, "Gluons", "PLE");
    leg1->AddEntry(quark_pt_resp, "Quarks", "PLE");
    leg1->AddEntry(gluon_pt_resp_nGenJetPF_w, "Gluons weighted", "PLE");
    leg1->AddEntry(quark_pt_resp_nGenJetPF_w, "Quarks weighted", "PLE");
    leg1->AddEntry(gluon_pt_resp_nGenJetPFSig_w, "Gluons weighted (sig)", "PLE");
    leg1->AddEntry(quark_pt_resp_nGenJetPFSig_w, "Quarks weighted (sig)", "PLE");


    TLatex *tex = new TLatex(); tex->SetNDC();
    tex->SetTextSize(0.05); tex->SetTextColor(kBlack);
    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //Response gaus
    TCanvas* c7 = tdrCanvas("c7", h, 4, 11, kSquare);
    h->GetXaxis()->SetNoExponent();
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(gluon_pt_resp_gaus, "P", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPF_w_gaus, "P", kFullTriangleUp, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPFSig_w_gaus, "P", kFullTriangleDown, kOrange +7);
    tdrDraw(quark_pt_resp_gaus, "P", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPF_w_gaus, "P", kFullTriangleUp, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPFSig_w_gaus, "P", kFullTriangleDown, kCyan +2);
    gluon_pt_resp_gaus->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPF_w_gaus->SetMarkerSize(1.2);
    quark_pt_resp_gaus->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPF_w_gaus->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPFSig_w_gaus->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPFSig_w_gaus->SetMarkerSize(1.2);

    TLegend *leg3 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg3->AddEntry(gluon_pt_resp_gaus, "Gluons", "PLE");
    leg3->AddEntry(quark_pt_resp_gaus, "Quarks", "PLE");
    leg3->AddEntry(gluon_pt_resp_nGenJetPF_w_gaus, "Gluons weighted", "PLE");
    leg3->AddEntry(quark_pt_resp_nGenJetPF_w_gaus, "Quarks weighted", "PLE");
    leg3->AddEntry(gluon_pt_resp_nGenJetPFSig_w_gaus, "Gluons weighted (sig)", "PLE");
    leg3->AddEntry(quark_pt_resp_nGenJetPFSig_w_gaus, "Quarks weighted (sig)", "PLE");
 
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

/*
    //genJetMass
    TCanvas* c5 = new TCanvas("c5","c5",1400,700);
    c5->Divide(2,1);

    c5->cd(1);
    gluon_genJetMass->GetXaxis()->SetTitle("p_{T} (GeV)");
    gluon_genJetMass->GetYaxis()->SetTitle("genJetMass");
    gluon_genJetMass->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    gluon_genJetMass->Draw("lego2");

    c5->cd(2);
    quark_genJetMass->GetXaxis()->SetTitle("p_{T} (GeV");
    quark_genJetMass->GetYaxis()->SetTitle("p_{T} (GeV)");
    quark_genJetMass->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    quark_genJetMass->Draw("lego2");
    */
}
