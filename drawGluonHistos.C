#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>

void drawGluonHistos() {
    TFile* f = TFile::Open("outputGluonHistos.root");

    // Gluons
   TProfile* gluon_pt_resp;
   TProfile* gluon_pt_resp_nGenJetPF;
   TProfile* gluon_nGenPF;

   // Quarks
   TProfile* quark_pt_resp;
   TProfile* quark_pt_resp_nGenJetPF;
   TProfile* quark_nGenPF;

    f->GetObject("gluon_pt_resp", gluon_pt_resp);
    f->GetObject("gluon_pt_resp_nGenJetPF", gluon_pt_resp_nGenJetPF);
    f->GetObject("gluon_nGenPF", gluon_nGenPF);
    f->GetObject("quark_pt_resp", quark_pt_resp);
    f->GetObject("quark_pt_resp_nGenJetPF", quark_pt_resp_nGenJetPF);
    f->GetObject("quark_nGenPF", quark_nGenPF);

    TCanvas* c = new TCanvas("c","",900,700);

    gluon_pt_resp->SetMarkerColor(kBlack);
    gluon_pt_resp->SetLineColor(kBlack);
    gluon_pt_resp->SetLineWidth(2);

    gluon_pt_resp_nGenJetPF->SetMarkerColor(kGray +1);
    gluon_pt_resp_nGenJetPF->SetLineColor(kGray +1);
    gluon_pt_resp_nGenJetPF->SetLineWidth(2);

    quark_pt_resp->SetMarkerColor(kRed);
    quark_pt_resp->SetLineColor(kRed);
    quark_pt_resp->SetLineWidth(2);

    quark_pt_resp_nGenJetPF->SetMarkerColor(kYellow);
    quark_pt_resp_nGenJetPF->SetLineColor(kYellow);
    quark_pt_resp_nGenJetPF->SetLineWidth(2);


    gPad->SetLogx();
    gluon_pt_resp->SetAxisRange(30, 3500, "X");
    gluon_pt_resp->SetAxisRange(0.97, 1.03, "Y");
    gluon_pt_resp->GetXaxis()->SetTitle("jetPt");
    gluon_pt_resp->GetXaxis()->SetMoreLogLabels();
    gluon_pt_resp->GetXaxis()->SetNoExponent();

    gluon_pt_resp->Draw();
    quark_pt_resp->Draw("same");
    gluon_pt_resp_nGenJetPF->Draw("same");
    quark_pt_resp_nGenJetPF->Draw("same");

    TCanvas* c2 = new TCanvas("c2", "", 900, 700);

    gluon_nGenPF->SetMarkerColor(kBlack);
    gluon_nGenPF->SetLineColor(kBlack);
    gluon_nGenPF->SetLineWidth(2);

    quark_nGenPF->SetMarkerColor(kRed);
    quark_nGenPF->SetLineColor(kRed);
    quark_nGenPF->SetLineWidth(2);

    gPad->SetLogx();
    gluon_nGenPF->SetAxisRange(30, 3500, "X");
    gluon_nGenPF->GetXaxis()->SetTitle("jetPt");
    gluon_nGenPF->GetXaxis()->SetMoreLogLabels();
    gluon_nGenPF->GetXaxis()->SetNoExponent();

    gluon_nGenPF->Draw();
    quark_nGenPF->Draw("same");
}
