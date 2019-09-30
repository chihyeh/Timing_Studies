#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array

f1= ROOT.TFile.Open("tev5mm_pythia6_zprime5tev_ww_1GeV_cut.root",'r')
f2= ROOT.TFile.Open("tev5mm_pythia6_zprime5tev_ww_1GeV_cut.root",'r')
h1 = f1.Get("Check_matching_0P2")
h2 = f2.Get("Check_matching_0P4")
h3 = f1.Get("Timing_Standard")

#h1 = f1.Get("Timing_detector_Trailing")
#h2 = f2.Get("Timing_detector_Trailing")

Data_Entries_1=h1.GetEntries()
Data_Entries_2=h2.GetEntries()
#h1.Sumw2()
#h1.Scale(1/h1.Integral())
#h2.Sumw2()
#h2.Scale(1/h2.Integral())
h3.Sumw2()
h3.Scale(1/h3.Integral())


a = TH1F ("a","a",10,0,100)
a.Fill(1)

c = TCanvas("c1", "c1",0,0,500,500)
gStyle.SetOptStat(0)

leg = TLegend(0.1,0.7,0.3,0.9)
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

h3.SetLineColor(3)
h3.SetLineWidth(3)
h3.SetLineStyle(1)

h1.SetMarkerStyle(9)
h2.SetMarkerStyle(9)
h3.SetMarkerStyle(9)

h2.GetXaxis().SetRangeUser(0,2)
h2.GetYaxis().SetRangeUser(0,4000)


h2.SetTitle("Jet number(s)")
h2.SetTitle("Jet number(s)")
h2.SetXTitle("Number of jet(s)")
h2.SetXTitle("Number of jet(s)")
h2.SetYTitle("Event number")
h2.SetYTitle("Event number")
leg.AddEntry("","FD group - SiFCC","")
#leg.AddEntry("","Z'(5TeV)#rightarrowq#bar{q}#rightarrow1 subjet","")
leg.AddEntry("","Z'(5TeV)#rightarrowW^{+}W^{-}#rightarrow2 subjet(s)","")
leg.AddEntry(h1,"Delta(R) < 0.2","l")
leg.AddEntry(h2,"Delta(R) < 0.4","l")

leg.Draw()

#Z'("+str(energy_array[1][m])+"TeV)#rightarrowt#bar{t}#rightarrow3 jet
#Z'("+str(energy_array[1][m])+"TeV)#rightarrowq#bar{q}#rightarrow1 jet
#Z'("+str(energy_array[1][m])+"TeV)#rightarrowW^{+}W^{-}#rightarrow2 jet
h1.GetYaxis().SetLabelSize(0.03)
h2.GetYaxis().SetLabelSize(0.03)
h1.GetXaxis().SetTitleFont(22)
h2.GetXaxis().SetTitleFont(22)
h1.GetYaxis().SetTitleFont(22)
h2.GetYaxis().SetTitleFont(22)
h1.GetXaxis().SetLabelFont(22)
h2.GetXaxis().SetLabelFont(22)
h1.GetYaxis().SetLabelFont(22)
h2.GetYaxis().SetLabelFont(22)

h2.Draw("hist")
h1.Draw("histsame")

leg.Draw()

c.Print("WW_checking_1GeV_cut.pdf")







