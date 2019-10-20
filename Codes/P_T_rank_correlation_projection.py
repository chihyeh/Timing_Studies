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
h1 = f1.Get("Timing_P_rank_difference_momentum_correlation")
h2 = f1.Get("Timing_P_rank_difference_momentum_correlation")


c = TCanvas("c1", "c1",0,0,500,500)
gStyle.SetOptStat(0)
cutg1 = TCutG("cutg1",5)
cutg1.SetPoint(0,0,-100)
cutg1.SetPoint(1,0,100)
cutg1.SetPoint(2,0.5,100)
cutg1.SetPoint(3,0.5,-100)
cutg1.SetPoint(4,0,-100)

cutg2 = TCutG("cutg2",5)
cutg2.SetPoint(0,0.5,-100)
cutg2.SetPoint(1,0.5,100)
cutg2.SetPoint(2,1,100)
cutg2.SetPoint(3,1,-100)
cutg2.SetPoint(4,0.5,-100)

cutg3 = TCutG("cutg3",5)
cutg3.SetPoint(0,1,-100)
cutg3.SetPoint(1,1,100)
cutg3.SetPoint(2,1.5,100)
cutg3.SetPoint(3,1.5,-100)
cutg3.SetPoint(4,1,-100)

cutg4 = TCutG("cutg4",5)
cutg4.SetPoint(0,1.5,-100)
cutg4.SetPoint(1,1.5,100)
cutg4.SetPoint(2,2,100)
cutg4.SetPoint(3,2,-100)
cutg4.SetPoint(4,1.5,-100)

cutg5 = TCutG("cutg5",5)
cutg5.SetPoint(0,2,-100)
cutg5.SetPoint(1,2,100)
cutg5.SetPoint(2,2.5,100)
cutg5.SetPoint(3,2.5,-100)
cutg5.SetPoint(4,2,-100)

cutg6 = TCutG("cutg6",5)
cutg6.SetPoint(0,2.5,-100)
cutg6.SetPoint(1,2.5,100)
cutg6.SetPoint(2,3,100)
cutg6.SetPoint(3,3,-100)
cutg6.SetPoint(4,2.5,-100)

H1_Projection = h1.ProjectionY("1",1,100,"[cutg1]")
H2_Projection = h1.ProjectionY("2",1,100,"[cutg2]")
H3_Projection = h1.ProjectionY("3",1,100,"[cutg3]")
H4_Projection = h1.ProjectionY("4",1,100,"[cutg4]")
H5_Projection = h1.ProjectionY("5",1,100,"[cutg5]")
H6_Projection = h1.ProjectionY("6",1,100,"[cutg6]")

H1_Projection.Sumw2()
H1_Projection.Scale(1/H1_Projection.Integral())

H2_Projection.Sumw2()
H2_Projection.Scale(1/H2_Projection.Integral())

H3_Projection.Sumw2()
H3_Projection.Scale(1/H3_Projection.Integral())

H4_Projection.Sumw2()
H4_Projection.Scale(1/H4_Projection.Integral())

H5_Projection.Sumw2()
H5_Projection.Scale(1/H5_Projection.Integral())

H6_Projection.Sumw2()
H6_Projection.Scale(1/H6_Projection.Integral())


H1_Projection.SetLineColor(1)
H1_Projection.SetLineWidth(2)
H1_Projection.SetLineStyle(1)

H2_Projection.SetLineColor(2)
H2_Projection.SetLineWidth(2)
H2_Projection.SetLineStyle(1)

H3_Projection.SetLineColor(3)
H3_Projection.SetLineWidth(2)
H3_Projection.SetLineStyle(1)

H4_Projection.SetLineColor(4)
H4_Projection.SetLineWidth(2)
H4_Projection.SetLineStyle(1)

H5_Projection.SetLineColor(5)
H5_Projection.SetLineWidth(2)
H5_Projection.SetLineStyle(1)

H6_Projection.SetLineColor(6)
H6_Projection.SetLineWidth(2)
H6_Projection.SetLineStyle(1)

leg = TLegend(0.1,0.6,0.4,0.9)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.SetBorderSize(0)
leg.SetTextFont(22)
leg.AddEntry("","FD group - SiFCC","")
leg.AddEntry("","Log(P) range:","")
leg.AddEntry(H1_Projection,"0-0.5","l")
leg.AddEntry(H2_Projection,"0.5-1","l")
leg.AddEntry(H3_Projection,"1-1.5","l")
leg.AddEntry(H4_Projection,"1.5-2","l")
leg.AddEntry(H5_Projection,"2-2.5","l")
leg.AddEntry(H6_Projection,"2.5-3","l")


H1_Projection.GetXaxis().SetRangeUser(-50,50)
H1_Projection.GetYaxis().SetRangeUser(0,0.3)

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
H3_Projection.Draw("histsame")
H4_Projection.Draw("histsame")
H5_Projection.Draw("histsame")
H6_Projection.Draw("histsame")

leg.Draw()

c.Print("P_T_rank_correlation_QQ_plot_slice1.pdf")






