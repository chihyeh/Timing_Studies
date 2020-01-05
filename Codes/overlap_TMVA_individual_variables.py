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

Energies=array("i",[5,10,20,40])
for i in range(0,4):
    c = TCanvas("c1", "c1",0,0,500,500)
    f4 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Codes/TMVA_for_timing_dR_PT_"+str(Energies[i])+"TeV_reco.root",'r')
    f5 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Codes/TMVA_for_timing_dR_PT_T_"+str(Energies[i])+"TeV_reco.root",'r')
    f6 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/Codes/TMVA_for_timing_dR_PT_"+str(Energies[i])+"TeV_reco_track.root",'r')

    #=====================#Get the histogram from TDirectoryFile
    h1_1 = f4.Get("dataset")
    h2_1 = h1_1.Get("Method_BDT")
    h3_1 = h2_1.Get("BDT")
    h4_1 = h3_1.Get("MVA_BDT_rejBvsS")
    h1_2 = f5.Get("dataset")
    h2_2 = h1_2.Get("Method_BDT")
    h3_2 = h2_2.Get("BDT")
    h4_2 = h3_2.Get("MVA_BDT_rejBvsS")
    h1_3 = f6.Get("dataset")
    h2_3 = h1_3.Get("Method_BDT")
    h3_3 = h2_3.Get("BDT")
    h4_3 = h3_3.Get("MVA_BDT_rejBvsS")

    #======================#

    h4_1.SetLineColor(1)
    h4_2.SetLineColor(7)
    h4_3.SetLineColor(2)

    gStyle.SetOptStat(0)
    h4_1.GetXaxis().SetRangeUser(0,1)
    h4_1.GetYaxis().SetRangeUser(0,1.1)
    h4_1.SetYTitle("1-Background efficiency")

    leg = TLegend(0.1,0.3,0.4,0.5)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetTextFont(22)
    leg.AddEntry("",str(Energies[i])+"TeV","")
    leg.AddEntry(h4_1,"TMVA(BDT)_dR_PT_ECAL","l")
    leg.AddEntry(h4_2,"TMVA(BDT)_dR_PT_T_ECAL","l")
    leg.AddEntry(h4_3,"TMVA(BDT)_dR_PT_track","l")

    h4_1.Draw("L")
    h4_2.Draw("Lsame")
    h4_3.Draw("Lsame")

    leg.Draw()
    c.Draw()

    c.Print("/Users/ms08962476/singularity/TIming_Studies/Codes/"+str(Energies[i])+"TeV_Reco/BDT_plot_dR_dRplusID_"+str(Energies[i])+"TeV_Reco.pdf")
