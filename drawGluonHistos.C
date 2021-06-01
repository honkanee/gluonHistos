#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>

void drawGluonHistos() {
    TFile* f = TFile::Open("outputGluonHistos.root");

    TProfile* jetPt_genJetPt_ratio;
    TProfile* jetPt_genJetPt_ratio_w_PF;

    f->GetObject("jetPt_genJetPt", jetPt_genJetPt_ratio);
    f->GetObject("jetPt_genJetPt, nGenJetPF", jetPt_genJetPt_ratio_w_PF);

    jetPt_genJetPt_ratio->SetMarkerColor(kBlack);
    jetPt_genJetPt_ratio->SetLineColor(kBlack);
    jetPt_genJetPt_ratio->SetLineWidth(2);

    jetPt_genJetPt_ratio_w_PF->SetMarkerColor(kRed);
    jetPt_genJetPt_ratio_w_PF->SetLineColor(kRed);
    jetPt_genJetPt_ratio_w_PF->SetLineWidth(2);

    TCanvas* c = new TCanvas("c","",900,700);

    jetPt_genJetPt_ratio->SetAxisRange(0.97, 1.03, "Y");
    jetPt_genJetPt_ratio->GetXaxis()->SetTitle("jetPt");

    jetPt_genJetPt_ratio->Draw();
    jetPt_genJetPt_ratio_w_PF->Draw("same");
}
