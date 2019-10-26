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



#Set eta==0 ~ eta==1
#Set momentum==0~20GeV
#Consider the effect of the magnetic field
#Resolution = 10,20,30,100,500,1000
def Theta_transformation(eta):
    return 2 * (math.atan(math.pow(2.718,-eta)))

def Hypothesis_particle_energy(Mass,momentum):
    return math.pow(math.pow(momentum,2)+math.pow(Mass,2),0.5)

def Resolution_term(momentum,theta,Resolution_time):
    Resolution_time_pico = Resolution_time * math.pow(10,-12)
    Lz_length = ((2.3) / (math.tan(theta)))
    P_Z = momentum * math.cos(theta)
    C = 3 * (math.pow(10,8))
    return (3*Resolution_time_pico*C*P_Z)/(Lz_length)

def mass_threshold(A,B,momentum):
    return math.pow(  (math.pow((A)+(B),2) - math.pow(momentum,2)) ,0.5)


yarray=array("f",[])
xarray=array("f",[])
xarray1=array("f",[])
xarray2=array("f",[])
xarray3=array("f",[])
xarray4=array("f",[])
xarray5=array("f",[])

xarrayerror=array("f",[])
yarrayerror=array("f",[])

#h1 = TF1("f1","[0]*x*sin([1]*x)",-3,3);
for i in range(1,101):
    Pion_mass = 0.134 #GeV
    Kaon_mass = 0.493 #GeV
    #Resolution = 10,20,30,100,500,1000
    Particle_momentum = (0.25*i)
    yarray.append(math.log10(Particle_momentum))
    #==============================================================================================================#
    M  =mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),10), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray.append(M)
    M_1=mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),20), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray1.append(M_1)
    M_2=mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),30), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray2.append(M_2)
    M_3=mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),100), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray3.append(M_3)
    M_4=mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),500), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray4.append(M_4)
    M_5=mass_threshold( Resolution_term(Particle_momentum,Theta_transformation(0),1000), Hypothesis_particle_energy(Kaon_mass,Particle_momentum),Particle_momentum )
    xarray5.append(M_5)
    #==============================================================================================================#
    xarrayerror.append(0)
    yarrayerror.append(0)


c = TCanvas("c1", "c1",0,0,500,500)
gStyle.SetOptStat(0)


gr  = TGraphErrors(100,xarray, yarray,xarrayerror,yarrayerror)
gr1 = TGraphErrors(100,xarray1,yarray,xarrayerror,yarrayerror)
gr2 = TGraphErrors(100,xarray2,yarray,xarrayerror,yarrayerror)
gr3 = TGraphErrors(100,xarray3,yarray,xarrayerror,yarrayerror)
gr4 = TGraphErrors(100,xarray4,yarray,xarrayerror,yarrayerror)
gr5 = TGraphErrors(100,xarray5,yarray,xarrayerror,yarrayerror)

gr.SetLineStyle(1)
gr.SetMarkerColor(2)
gr.SetMarkerStyle(8)
gr.SetLineWidth(2)

gr1.SetLineStyle(1)
gr1.SetMarkerColor(3)
gr1.SetMarkerStyle(8)
gr1.SetLineWidth(2)

gr2.SetLineStyle(1)
gr2.SetMarkerColor(4)
gr2.SetMarkerStyle(8)
gr2.SetLineWidth(2)

gr3.SetLineStyle(1)
gr3.SetMarkerColor(5)
gr3.SetMarkerStyle(8)
gr3.SetLineWidth(2)

gr4.SetLineStyle(1)
gr4.SetMarkerColor(6)
gr4.SetMarkerStyle(8)
gr4.SetLineWidth(2)

gr5.SetLineStyle(1)
gr5.SetMarkerColor(7)
gr5.SetMarkerStyle(8)
gr5.SetLineWidth(2)

gr.SetTitle(";mass[GeV] ;Log_{10}P  [GeV] ")
gr.GetXaxis().SetRangeUser(0,10)
gr.GetYaxis().SetRangeUser(0,2)
#gr.GetHistogram().SetMaximum(6)
#gr.GetHistogram().SetMinimum(-2)
#gr.GetXaxis().SetLimits(0,1000)
gr.GetXaxis().CenterTitle()
gr.GetYaxis().CenterTitle()
gr.GetXaxis().SetTitleOffset(1)
gr.GetYaxis().SetTitleOffset(1)

#================================
gr.GetYaxis().SetTitleSize(0.04)
gr.GetXaxis().SetTitleSize(0.04)
gr.GetXaxis().SetLabelSize(0.04)
gr.GetYaxis().SetLabelSize(0.04)
gr.GetXaxis().SetLabelFont(22)
gr.GetYaxis().SetLabelFont(22)
gr.GetXaxis().SetTitleColor(1)
gr.GetYaxis().SetTitleColor(1)

leg1=TLegend(0.2,0.6,0.4,0.9)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.05)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry("","FD group - SiFCC","")
#leg.AddEntry(h1,"Z'(5TeV)#rightarrowq#bar{q}#rightarrow1 subjet","l")
#leg.AddEntry(h2,"Z'(5TeV)#rightarrowq#bar{q}#rightarrow1 subjet(#eta cut)","l")
leg1.AddEntry(gr, "10ps","lp")
leg1.AddEntry(gr1,"20ps","lp")
leg1.AddEntry(gr2,"30ps","lp")
leg1.AddEntry(gr3,"100ps","lp")
leg1.AddEntry(gr4,"500ps","lp")
leg1.AddEntry(gr5,"1000ps","lp")

c.Draw()
gr.Draw("ALP")
gr1.Draw("LPsame")
gr2.Draw("LPsame")
gr3.Draw("LPsame")
gr4.Draw("LPsame")
gr5.Draw("LPsame")
leg1.Draw()
c.Print("3_sigma_kaon_hypothesis_eta_0.pdf")










