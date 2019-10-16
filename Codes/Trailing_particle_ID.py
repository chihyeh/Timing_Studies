#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors,TText
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array

f1= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Files/tev5mm_pythia6_zprime5tev_qq.root",'r')
f2= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Files/tev5mm_pythia6_zprime5tev_ww.root",'r')
f3= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_qq_1P5GeV_cut.root",'r')
f4= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_ww_1P5GeV_cut.root",'r')

h1 = f1.Get("Timing_detector_next_to_trailing")
h2 = f2.Get("Timing_detector_next_to_trailing")
h3 = f3.Get("Trailing_particle_ID")
h4 = f4.Get("Trailing_particle_ID")
h5 = f1.Get("Timing_Standard")

h1.Sumw2()
h1.Scale(1/h1.Integral())
h2.Sumw2()
h2.Scale(1/h2.Integral())
h3.Sumw2()
h3.Scale(1/h3.Integral())
h4.Sumw2()
h4.Scale(1/h4.Integral())
h5.Sumw2()
h5.Scale(1/h5.Integral())


a = TH1F ("a","a",10,0,100)
a.Fill(-1)

c = TCanvas("c1", "c1",0,0,500,500)
gStyle.SetOptStat(0)

leg = TLegend(0.35,0.75,0.65,0.9)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.SetBorderSize(0)
leg.SetTextFont(22)
leg.Draw()

h1.SetLineColor(1)
h1.SetLineWidth(2)
h1.SetLineStyle(1)

h2.SetLineColor(2)
h2.SetLineWidth(2)
h2.SetLineStyle(1)

h3.SetLineColor(1)
h3.SetLineWidth(2)
h3.SetLineStyle(1)

h4.SetLineColor(2)
h4.SetLineWidth(2)
h4.SetLineStyle(1)

h5.SetLineColor(3)
h5.SetLineWidth(2)
h5.SetLineStyle(1)

h1.SetMarkerStyle(9)
h2.SetMarkerStyle(9)
h3.SetMarkerStyle(9)
h4.SetMarkerStyle(9)
h5.SetMarkerStyle(9)

h3.GetXaxis().SetRangeUser(0,16)
h3.GetYaxis().SetRangeUser(0,0.6)
h3.GetYaxis().SetRangeUser(0,0.6)


t =  TLatex(0.3,.05,"e^{-}");
t1 =  TLatex(1.3,.05,"#nu_{e}");
t2 =  TLatex(2.3,.05,"#mu^{-}");
t3 =  TLatex(3.3,.05,"#nu_{#mu}");
t4 =  TLatex(4.3,.05,"#gamma");
t5 =  TLatex(5.1,.05,"K_{L}^{0}");
t6 =  TLatex(6.3,.05,"#pi^{+}");
t7 =  TLatex(7.1,.05,"K_{S}^{0}");
t8 =  TLatex(8.1,.05,"K^{+}");
t9 =  TLatex(9.3,.05,"n");
t10 =  TLatex(10.3,.05,"p");
t11 =  TLatex(11.3,.05,"#Sigma^{-}");
t12 =  TLatex(12.3,.05,"#Lambda");
t13 =  TLatex(13.3,.05,"#Xi^{-}");
t14 =  TLatex(14.3,.05,"#Sigma^{+}");
t15 =  TLatex(15.3,.05,"#Xi^{0}");

h3.SetTitle("Trailing particle ID")
h3.SetTitle("Trailing particle ID")
h3.SetXTitle("Kinds of particles")
h3.SetXTitle("Kinds of particles")
h3.SetYTitle("Arbitrary number")
h3.SetYTitle("Arbitrary number")
leg.AddEntry("","FD group - SiFCC","")
leg.AddEntry(h3,"Z'(5TeV)#rightarrowq#bar{q}#rightarrow1 subjet","l")
leg.AddEntry(h4,"Z'(5TeV)#rightarrowW^{+}W^{-}#rightarrow2 subjets","l")

leg.Draw()

#Z'("+str(energy_array[1][m])+"TeV)#rightarrowt#bar{t}#rightarrow3 jet
#Z'("+str(energy_array[1][m])+"TeV)#rightarrowq#bar{q}#rightarrow1 jet
#Z'("+str(energy_array[1][m])+"TeV)#rightarrowW^{+}W^{-}#rightarrow2 jet
h3.GetYaxis().SetLabelSize(0.03)
h4.GetYaxis().SetLabelSize(0.03)
h3.GetXaxis().SetTitleFont(22)
h4.GetXaxis().SetTitleFont(22)
h3.GetYaxis().SetTitleFont(22)
h4.GetYaxis().SetTitleFont(22)
h3.GetXaxis().SetLabelFont(22)
h4.GetXaxis().SetLabelFont(22)
h3.GetYaxis().SetLabelFont(22)
h4.GetYaxis().SetLabelFont(22)

h3.Draw("hist")
h4.Draw("histsame")
t.Draw("same")
t1.Draw("same")
t2.Draw("same")
t3.Draw("same")
t4.Draw("same")
t5.Draw("same")
t6.Draw("same")
t7.Draw("same")
t8.Draw("same")
t9.Draw("same")
t10.Draw("same")
t11.Draw("same")
t12.Draw("same")
t13.Draw("same")
t14.Draw("same")
t15.Draw("same")

leg.Draw()

c.Print("Try_Trailing_particles_ID.pdf")







