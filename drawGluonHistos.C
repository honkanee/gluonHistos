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

    TH2D* gluon_nGenPF_probs;
    TH2D* gluon_genJetMass;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF_prof", gluon_nGenPF);
    f->GetObject("gluon_nGenPF_prof_w", gluon_nGenPF_w);
    f->GetObject("gluon_pt_resp_nGenJetPF_w", gluon_pt_resp_nGenJetPF_w);
    f->GetObject("gluon_genJetMass", gluon_genJetMass);

    f2->GetObject("gluon_probs_smooth", gluon_nGenPF_probs);

    // Quarks
    TProfile* quark_pt_resp;
    TProfile* quark_pt_resp_nGenJetPF;
    TProfile* quark_nGenPF;
    TProfile* quark_nGenPF_w;
    TProfile* quark_pt_resp_nGenJetPF_w;

    TH2D* quark_nGenPF_probs;
    TH2D* quark_genJetMass;
    TH2D* nUEPF_probs;

    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF_prof", quark_nGenPF);
    f->GetObject("quark_nGenPF_prof_w", quark_nGenPF_w);
    f->GetObject("quark_pt_resp_nGenJetPF_w", quark_pt_resp_nGenJetPF_w);
    f->GetObject("quark_genJetMass", quark_genJetMass);

    f2->GetObject("quark_probs_smooth", quark_nGenPF_probs);

    TH2D* nGenJetPF_likelyhood;

    f->GetObject("nGenJetPF_likelyhood", nGenJetPF_likelyhood);
    f2->GetObject("UE_nPF_probs", nUEPF_probs);

    //Pt response
    TH1D* h = tdrHist("h", "Response", 0.9, 1.2,
    "p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c = tdrCanvas("c", h, 4, 11, kSquare);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(gluon_pt_resp, "", kFullCircle, kRed);
    tdrDraw(gluon_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kRed);
    tdrDraw(quark_pt_resp, "", kFullCircle, kBlue);
    tdrDraw(quark_pt_resp_nGenJetPF_w, "", kFullTriangleUp, kBlue);
    gluon_pt_resp->SetMarkerSize(1.2);
    gluon_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);
    quark_pt_resp->SetMarkerSize(1.2);
    quark_pt_resp_nGenJetPF_w->SetMarkerSize(1.2);

    TLegend *leg1 = tdrLeg(0.37,0.90-6*0.045,0.57,0.90);
    leg1->AddEntry(gluon_pt_resp, "Gluons", "PLE");
    leg1->AddEntry(quark_pt_resp, "Quarks", "PLE");
    leg1->AddEntry(gluon_pt_resp_nGenJetPF_w, "Gluons weighted", "PLE");
    leg1->AddEntry(quark_pt_resp_nGenJetPF_w, "Quarks weighted", "PLE");

    TLatex *tex = new TLatex(); tex->SetNDC();
    tex->SetTextSize(0.05); tex->SetTextColor(kBlack);
    tex->DrawLatex(0.15, 0.95,"|#eta|<1.3");

    //nGenJetPF profiles
    TH1D* h2 = tdrHist("h2", "nGenPF", 10, 70,
    "p_{T} (GeV)", 30, 3500);
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
    gluon_nGenPF_probs->GetXaxis()->SetTitle("pT (GeV)");
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
    quark_nGenPF_probs->GetXaxis()->SetTitle("p_{T} (GeV)");
    quark_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPF_probs->GetZaxis()->SetTitle("P");
    quark_nGenPF_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPF_probs->SetTitle("quarks");
    quark_nGenPF_probs->Draw("lego2");

    
    TText *qt = new TText(0.5, 0.5, "quarks");
    qt->Draw();

    //nGenJet probablity distributions 2D
    TCanvas* c4 = new TCanvas("c4","c4",1400,700);
    c4->Divide(2,1);

    c4->cd(1);
    gluon_nGenPF_probs->SetAxisRange(30, 3500, "X");
    gPad->SetLogx();
    gluon_nGenPF_probs->GetXaxis()->SetTitle("pT (GeV)");
    gluon_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    gluon_nGenPF_probs->GetZaxis()->SetTitle("prob");
    gluon_nGenPF_probs->SetAxisRange(0,0.1,"Z");
    gluon_nGenPF_probs->Draw("colz2");

    gt->Draw();

    c4->cd(2);
    quark_nGenPF_probs->SetAxisRange(30, 3500,"X");
    gPad->SetLogx();
    quark_nGenPF_probs->GetXaxis()->SetTitle("pT (GeV)");
    quark_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPF_probs->GetZaxis()->SetTitle("prob");
    quark_nGenPF_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPF_probs->Draw("colz2");

    qt->Draw();

    TH1D* h3 = tdrHist("h3", "nGenPF", 10, 70,
    "p_{T} (GeV)", 30, 3500);
    h3->GetXaxis()->SetNoExponent();
    h3->GetXaxis()->SetMoreLogLabels();
    tdrCanvas("c6", h3, 4, 11,  kSquare);
    gPad->SetLogx();

    nGenJetPF_likelyhood->SetAxisRange(30, 3500, "X");
    nGenJetPF_likelyhood->SetAxisRange(0, 150, "Y");
    nGenJetPF_likelyhood->GetXaxis()->SetTitle("p_{T} (GeV)");
    nGenJetPF_likelyhood->GetYaxis()->SetTitle("n_{GenJetPF}");
    nGenJetPF_likelyhood->Draw("colz");

    TH1D* h4 = tdrHist("h4", "nUEPF", 0, 7,
    "p_{T} (GeV)", 30, 3500);
    tdrCanvas("c7", h4, 4, 11,  kSquare);
    gPad->SetLogx();
    h4->GetXaxis()->SetNoExponent();
    h4->GetXaxis()->SetMoreLogLabels();

    nUEPF_probs->SetAxisRange(30, 3500, "X");
    nUEPF_probs->GetXaxis()->SetTitle("pT");
    nUEPF_probs->GetYaxis()->SetTitle("nUEPF");
    nUEPF_probs->GetZaxis()->SetTitle("P");
    nUEPF_probs->SetAxisRange(0, 150, "Y");
    nUEPF_probs->Draw("lego2");


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
