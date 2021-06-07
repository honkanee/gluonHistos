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

    // Gluons
    TProfile* gluon_pt_resp;
    TProfile* gluon_pt_resp_nGenJetPF;
    TProfile* gluon_nGenPF;
    TProfile* gluon_nGenPF_w;

    TH2D* gluon_nGenPF_probs;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF_prof", gluon_nGenPF);
    f->GetObject("gluon_nGenPF_prof_w", gluon_nGenPF_w);
    f->GetObject("gluon_nGenPF_probs", gluon_nGenPF_probs);

    // Quarks
    TProfile* quark_pt_resp;
    TProfile* quark_pt_resp_nGenJetPF;
    TProfile* quark_nGenPF;
    TProfile* quark_nGenPF_w;

    TH2D* quark_nGenPF_probs;

    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF_prof", quark_nGenPF);
    f->GetObject("quark_nGenPF_prof_w", quark_nGenPF_w);
    f->GetObject("quark_nGenPF_probs", quark_nGenPF_probs);

    //Pt response
    TCanvas* c = new TCanvas("c","",900,700);

    gluon_pt_resp->SetMarkerColor(kBlack);
    gluon_pt_resp->SetLineColor(kBlack);
    gluon_pt_resp->SetLineWidth(2);

    quark_pt_resp->SetMarkerColor(kRed);
    quark_pt_resp->SetLineColor(kRed);
    quark_pt_resp->SetLineWidth(2);

    gPad->SetLogx();
    gluon_pt_resp->SetAxisRange(30, 3500, "X");
    gluon_pt_resp->SetAxisRange(0.97, 1.03, "Y");
    gluon_pt_resp->GetXaxis()->SetTitle("jetPt");
    gluon_pt_resp->GetXaxis()->SetMoreLogLabels();
    gluon_pt_resp->GetXaxis()->SetNoExponent();

    gluon_pt_resp->Draw();
    quark_pt_resp->Draw("same");


    //nGenJetPF profiles
    TH1D* h = tdrHist("h", "nGenPF", 10, 70,
    "p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c2 = tdrCanvas("c2", h, 4);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    tdrDraw(gluon_nGenPF,"same", kFullCircle, kRed);
    tdrDraw(quark_nGenPF, "same", kOpenCircle, kRed -5);
    tdrDraw(gluon_nGenPF_w, "same", kFullTriangleUp, kBlue);
    tdrDraw(quark_nGenPF_w, "same", kOpenTriangleUp, kBlue -5);

    TLegend *leg1 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(gluon_nGenPF, "Gluons nGenPF");
    leg1->AddEntry(quark_nGenPF, "Quarks nGenPF");
    leg1->AddEntry(gluon_nGenPF_w, "Gluons weighted");
    leg1->AddEntry(quark_nGenPF_w, "Quarks weighted");

    TLatex *tex = new TLatex(); tex->SetNDC();
    tex->SetTextSize(0.04); tex->SetTextColor(kBlack);
    tex->DrawLatex(0.15, 0.6,"|#eta|<1.3");

    //Probablity distributions
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
    quark_nGenPF_probs->GetXaxis()->SetTitle("pT (GeV)");
    quark_nGenPF_probs->GetYaxis()->SetTitle("nGenPF");
    quark_nGenPF_probs->GetZaxis()->SetTitle("prob");
    quark_nGenPF_probs->SetAxisRange(0,0.1, "Z");
    quark_nGenPF_probs->SetTitle("quarks");
    quark_nGenPF_probs->Draw("lego2");

    
    TText *qt = new TText(0.5, 0.5, "quarks");
    qt->Draw();
}
