#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors
from ROOT import TH1D, TH1, TH1I, TCut, TCutG
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array

f1= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_qq_1P5GeV_cut_rank_reduce_tosix_mass.root",'r')
f2= ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_ww_1P5GeV_cut_rank_reduce_tosix_mass.root",'r')

h1 = f1.Get("Timing_P_rank_difference_momentum_correlation")
h2 = f2.Get("Timing_P_rank_difference_momentum_correlation")

for i in range(6):
    c = TCanvas("c1", "c1",0,0,500,500)
    gStyle.SetOptStat(0)
    cutg1 = TCutG("cutg1",5)
    cutg1.SetPoint(0,0.5*(i),-100)
    cutg1.SetPoint(1,0.5*(i),100)
    cutg1.SetPoint(2,0.5*(1+i),100)
    cutg1.SetPoint(3,0.5*(1+i),-100)
    cutg1.SetPoint(4,0.5*(i),-100)


    H1_Projection = h1.ProjectionY("1",1,100,"[cutg1]")
    H2_Projection = h2.ProjectionY("2",1,100,"[cutg1]")

    H1_Projection.Sumw2()
    H1_Projection.Scale(1/H1_Projection.Integral())

    H2_Projection.Sumw2()
    H2_Projection.Scale(1/H2_Projection.Integral())

    H1_Projection.SetLineColor(1)
    H1_Projection.SetLineWidth(2)
    H1_Projection.SetLineStyle(1)

    H2_Projection.SetLineColor(2)
    H2_Projection.SetLineWidth(2)
    H2_Projection.SetLineStyle(1)

    leg = TLegend(0.1,0.7,0.4,0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetTextFont(22)
    leg.AddEntry("","FD group - SiFCC","")
    leg.AddEntry(H1_Projection,"Z'(5TeV)#rightarrowq#bar{q}#rightarrow1 subjet","l")
    leg.AddEntry(H2_Projection,"Z'(5TeV)#rightarrowW^{+}W^{-}#rightarrow2 subjets","l")


    H1_Projection.GetXaxis().SetRangeUser(-50,50)
    H1_Projection.GetYaxis().SetRangeUser(0,0.3)
    H1_Projection.SetTitle("Log(P): "+str(0.5*(i))+"~"+str(0.5*(i+1))+" [GeV]")
    H1_Projection.SetXTitle("T rank-P rank")
    H1_Projection.SetYTitle("Arbitrary number")

    H1_Projection.GetXaxis().SetTitleOffset(1)
    H1_Projection.GetYaxis().SetLabelSize(0.03)
    H1_Projection.GetXaxis().SetTitleFont(22)
    H1_Projection.GetYaxis().SetTitleFont(22)
    H1_Projection.GetXaxis().SetLabelFont(22)
    H1_Projection.GetYaxis().SetLabelFont(22)


    #h1.Draw("colz")
    H1_Projection.Draw("hist")
    H2_Projection.Draw("histsame")

    leg.Draw()

    c.Print("P_T_rank_correlation_QQ_plot_slice1_"+str(i)+".pdf")






