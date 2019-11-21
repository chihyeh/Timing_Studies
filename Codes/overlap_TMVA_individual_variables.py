import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors, TDirectoryFile
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array

c = TCanvas("c1", "c1",0,0,500,500)
f4 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Codes/TMVA_for_timing_dR_PT_5TeV_reco.root",'r')
f5 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Codes/TMVA_for_timing_dR_PT_T_5TeV_reco.root",'r')

#=====================#Get the histogram from TDirectoryFile
h1_1 = f4.Get("dataset")
h2_1 = h1_1.Get("Method_BDT")
h3_1 = h2_1.Get("BDT")
h4_1 = h3_1.Get("MVA_BDT_rejBvsS")
h1_2 = f5.Get("dataset")
h2_2 = h1_2.Get("Method_BDT")
h3_2 = h2_2.Get("BDT")
h4_2 = h3_2.Get("MVA_BDT_rejBvsS")

#======================#

h4_1.SetLineColor(1)
h4_2.SetLineColor(7)
gStyle.SetOptStat(0)
h4_1.GetXaxis().SetRangeUser(0,1)
h4_1.GetYaxis().SetRangeUser(0,1.1)


leg = TLegend(0.1,0.7,0.4,0.9)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.SetBorderSize(0)
leg.SetTextFont(22)

leg.AddEntry(h4_1,"TMVA(BDT)_dR_PT","l")
leg.AddEntry(h4_2,"TMVA(BDT)_dR_PT_T","l")

h4_1.Draw("L")
h4_2.Draw("Lsame")
leg.Draw()
c.Draw()

c.Print("/Users/ms08962476/singularity/TIming_Studies/Codes/5TeV_Reco/BDT_plot_dR_dRplusID_5TeV_Reco.pdf")
