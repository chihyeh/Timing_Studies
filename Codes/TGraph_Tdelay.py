import ROOT
import sys
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraphErrors, TGraphAsymmErrors,TMultiGraph,TText,TNamed, TLatex, TF1, TFormula
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array
f1 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_qq_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f2 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev10mm_pythia6_zprime10tev_qq_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f3 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev20mm_pythia6_zprime20tev_qq_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f4 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev40mm_pythia6_zprime40tev_qq_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f5 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev5mm_pythia6_zprime5tev_ww_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f6 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev10mm_pythia6_zprime10tev_ww_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f7 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev20mm_pythia6_zprime20tev_ww_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')
f8 = ROOT.TFile.Open("/Users/ms08962476/singularity/TIming_Studies/tev40mm_pythia6_zprime40tev_ww_with_Eta_cut_for_component_check_1_reco_Com.root", 'r')

c = TCanvas("c1", "c1",0,0,500,500)

f1F =f1.Get("Timing_detecto_ECAL_TDif")
f2F =f2.Get("Timing_detecto_ECAL_TDif")
f3F =f3.Get("Timing_detecto_ECAL_TDif")
f4F =f4.Get("Timing_detecto_ECAL_TDif")
f5F =f5.Get("Timing_detecto_ECAL_TDif")
f6F =f6.Get("Timing_detecto_ECAL_TDif")
f7F =f7.Get("Timing_detecto_ECAL_TDif")
f8F =f8.Get("Timing_detecto_ECAL_TDif")

xarray=array("f",[5,10,20,40])
yarray_QQ=array("f",[])
yarray_WW=array("f",[])
xarray_error=array("f",[0,0,0,0])
yarray_QQ_RMS=array("f",[])
yarray_WW_RMS=array("f",[])

#===============================#
yarray_QQ.append(f1F.GetMean())
yarray_QQ.append(f2F.GetMean())
yarray_QQ.append(f3F.GetMean())
yarray_QQ.append(f4F.GetMean())
yarray_WW.append(f5F.GetMean())
yarray_WW.append(f6F.GetMean())
yarray_WW.append(f7F.GetMean())
yarray_WW.append(f8F.GetMean())
#===============================#
yarray_QQ_RMS.append(f1F.GetMeanError())
yarray_QQ_RMS.append(f2F.GetMeanError())
yarray_QQ_RMS.append(f3F.GetMeanError())
yarray_QQ_RMS.append(f4F.GetMeanError())
yarray_WW_RMS.append(f5F.GetMeanError())
yarray_WW_RMS.append(f6F.GetMeanError())
yarray_WW_RMS.append(f7F.GetMeanError())
yarray_WW_RMS.append(f8F.GetMeanError())
#===============================#


c = TCanvas("c1", "c1",0,0,500,500)

gr_QQ = TGraphErrors(4,xarray,yarray_QQ,xarray_error,yarray_QQ_RMS)
gr_QQ.SetLineColor(1)
gr_QQ.SetLineWidth(1)
gr_QQ.SetLineStyle(1)
gr_QQ.SetMarkerColor(2)
gr_QQ.SetMarkerStyle(1)
gr_QQ.SetMarkerSize(1)
gr_QQ.GetXaxis().SetTitle("E[TeV]")
gr_QQ.GetYaxis().SetTitle("<T_{delay}>[ns]")

gr_WW = TGraphErrors(4,xarray,yarray_WW,xarray_error,yarray_WW_RMS)
gr_WW.SetLineColor(3)
gr_WW.SetLineWidth(1)
gr_WW.SetLineStyle(1)
gr_WW.SetMarkerColor(4)
gr_WW.SetMarkerStyle(1)
gr_WW.SetMarkerSize(1)
#gr.GetXaxis().SetTitle("E[TeV]")
#gr.GetXaxis().SetTitleColor(4)
#gr.GetYaxis().SetTitle("<T_{delay>}")

gr_QQ.SetTitle("<T_{delay}>[ns] Vs E[TeV]")
gr_QQ.Draw("ALP")
gr_WW.Draw("LPsame")


#=================================
'''
fitFcn = TF1("fitFcn","sqrt([0]*[0]/(x)+[1]*[1])",20,300)
fitFcn.SetLineColor(2)
gr.Fit(fitFcn,"R")
print fitFcn.GetChisquare()
print fitFcn.GetNDF()
print fitFcn.GetChisquare()/fitFcn.GetNDF()

fitFcn1 = TF1("fitFcn1","sqrt([0]*[0]/(x)+[1]*[1])",20,300)
fitFcn1.SetLineColor(3)
gr1.Fit(fitFcn1,"R")
print fitFcn1.GetChisquare()
print fitFcn1.GetNDF()
print fitFcn1.GetChisquare()/fitFcn1.GetNDF()
'''

gr_QQ.GetHistogram().SetMaximum(+5)
gr_QQ.GetHistogram().SetMinimum(-5)
#gr.GetXaxis().SetLimits(0,320)
gr_QQ.GetXaxis().CenterTitle()
gr_QQ.GetYaxis().CenterTitle()

#=================================
leg1=TLegend(0.15,0.75,0.35,0.9)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.04)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry(gr_QQ,"Z'#rightarrowQQ","lp")
leg1.AddEntry(gr_WW,"Z'#rightarrowWW","lp")

c.Draw()
leg1.Draw()
c.Print("TGraph_TDelay.pdf")




