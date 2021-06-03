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
   TH1D* gluon_nGenPF_w;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF", gluon_nGenPF);
    f->GetObject("gluon_nGenPF_w", gluon_nGenPF_w);

   // Quarks
   TProfile* quark_pt_resp;
   TProfile* quark_pt_resp_nGenJetPF;
   TProfile* quark_nGenPF;
   TH1D* quark_nGenPF_w;

    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF", quark_nGenPF);
    f->GetObject("quark_nGenPF_w", quark_nGenPF_w);


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


    //nGenJetPF
    TH1D* h = tdrHist("h", "nGenPF", 10, 70,
    "p_{T} (GeV)", 30, 3500);
    h->GetXaxis()->SetNoExponent();
    TCanvas* c2 = tdrCanvas("c2", h, 4);
    h->GetXaxis()->SetMoreLogLabels();
    gPad->SetLogx();

    gluon_nGenPF->SetMarkerColor(kBlue);
    gluon_nGenPF->SetLineColor(kBlue);
    gluon_nGenPF->SetLineWidth(2);

    gluon_nGenPF_w->SetMarkerColor(kBlue-5);
    gluon_nGenPF_w->SetLineColor(kBlue-5);
    gluon_nGenPF_w->SetLineWidth(2);

    quark_nGenPF->SetMarkerColor(kRed);
    quark_nGenPF->SetLineColor(kRed);
    quark_nGenPF->SetLineWidth(2);

    quark_nGenPF_w->SetMarkerColor(kRed-5);
    quark_nGenPF_w->SetLineColor(kRed-5);
    quark_nGenPF_w->SetLineWidth(2);

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
}
