import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--bitfieldEcal", help="bitfield of dd4hep ECAL readout", type = str)
calo_init.parser.add_argument("--bitfieldHcal", help="bitfield of dd4hep HCAL readout", type = str)
calo_init.parse_args()
calo_init.print_config()

Thr = 0.
# 0.993327,1,0.417989,-6.35996e-06 (200ev of all energies except 20GeV) # 0.978295,1,0.479623,-5.39949e-06 (10,100,1000,10000 GeV a 600ev)
A = 0.978295
B = 0.479623
C = -5.39949e-06
# Linearity corr after minimisation
aLin = -9.19881e-01
bLin = 9.50605e-01
cLin = 1.00594e+00

ecalBitfield = "system:4,cryo:1,type:3,subtype:3,layer:8,eta:9,phi:10"
hcalBitfield = "system:4,layer:5,eta:9,phi:10"
hcalVolBitfield = "system:4,module:8,row:9,layer:5"

#if calo_init.args.bitfieldEcal:
#    ecalBitfield = calo_init.args.bitfieldEcal
#if calo_init.args.bitfieldHcal:
#    hcalBitfield = calo_init.args.bitfieldHcal

from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
import ROOT as r
import numpy as n
from draw_functions import draw_1histogram, draw_2histograms, draw_hist2d
from export_to_csv import *

r.gROOT.SetBatch(True)

Energy = float(calo_init.energy(calo_init.filenamesIn))
Et = float(Energy / r.TMath.CosH(0.))
sigma = r.TMath.Sqrt(r.TMath.Power(0.8*r.TMath.Sqrt(Energy)+0.02*(Energy),2))

etaBinsECal=40
etaBinsHCal=20
phiBinsECal=40
phiBinsHCal=20
etaMinECal=.2
phiMinECal=.2
etaMin=.25
phiMin=.25
 
h_etaphiEnergy1 = r.TH2D("h_etaphiEnergy1","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy2 = r.TH2D("h_etaphiEnergy2","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy3 = r.TH2D("h_etaphiEnergy3","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy4 = r.TH2D("h_etaphiEnergy4","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy5 = r.TH2D("h_etaphiEnergy5","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy6 = r.TH2D("h_etaphiEnergy6","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy7 = r.TH2D("h_etaphiEnergy7","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy8 = r.TH2D("h_etaphiEnergy8","", etaBinsECal,-etaMinECal,etaMinECal, phiBinsECal, -phiMinECal, phiMinECal)
h_etaphiEnergy9 = r.TH2D("h_etaphiEnergy9","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy10 = r.TH2D("h_etaphiEnergy10","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy11 = r.TH2D("h_etaphiEnergy11","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy12 = r.TH2D("h_etaphiEnergy12","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy13 = r.TH2D("h_etaphiEnergy13","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy14 = r.TH2D("h_etaphiEnergy14","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy15 = r.TH2D("h_etaphiEnergy15","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy16 = r.TH2D("h_etaphiEnergy16","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy17 = r.TH2D("h_etaphiEnergy17","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)
h_etaphiEnergy18 = r.TH2D("h_etaphiEnergy18","", etaBinsHCal,-etaMin,etaMin, phiBinsHCal, -phiMin, phiMin)

s_EPhi = r.THStack("s_EPhi","")

import glob

inputName = "/eos/experiment/fcc/hh/simulation/samples/v03/singlePart/pim/bFieldOn/eta0.36/"+str(int(Energy))+"GeV/ntup/positions/resegmentedHCal/output_100113638.root"
print inputName

chain_data = r.TChain("events") 
chain_data.Add( inputName )

#event = r.Event("event")  #object must be created before setting the branch address      
#rec_eta = r.std.vector(float)()
#rec_phi = r.std.vector(float)()
#rec_ene = r.std.vector(float)()
#rec_layer = r.std.vector(int)()
#ev_num = n.zeros(1, dtype=int) 
#chain_data.SetBranchAddress("event", event);
#chain_data.SetBranchAddress("rechit_eta", rec_eta)
#chain_data.SetBranchAddress("rechit_phi", rec_phi)
#chain_data.SetBranchAddress("rechit_energy", rec_ene)
#chain_data.SetBranchAddress("rechit_layer", rec_layer)


for i in range(1, 18):
    histName = "h_etaphiEnergy"+str(i)
    print histName
    command = "(TVector2::Phi_mpi_pi(rechit_phi-gen_phi[0])):rechit_eta-gen_eta[0]>>"+str(histName)
    cut = "rechit_layer == "+str(i)+" && ev_num==0"
    
    chain_data.Draw(command, cut, "col")

#
#for i in range(0, 1):
#    chain_data.GetEvent(i)
#
#    command = "rechit_eta:rechit_phi:rechit_energy>>h_etaphiEnergy"+str(i)
#    cut = "rechit_layer == "+str(i-1)+" && ev_num==0"
#    for hit in range(0, rec_eta.length):
#        h_etaphiEnergy1.Fill(rec_eta[hit], rec_phi[hit], rec_energy[hit])

print "Histograms saved in rootFiles/bFieldOn/official/cells_resegmentedHCal_calibrated_pi-"+str(int(Energy))+"GeV_eta036.root"
rootFileName = ("rootFiles/bFieldOn/official/cells_resegmentedHCal_wNoise_calibrated_pi-"+str(int(Energy))+"GeV_eta036.root")
f = r.TFile(rootFileName,"recreate")

cstack = r.TCanvas("cstack","cstack")
s_EPhi.Draw()

vecHistEnergy = []

vecHistEnergy.append(h_etaphiEnergy1)
vecHistEnergy.append(h_etaphiEnergy2)
vecHistEnergy.append(h_etaphiEnergy3)
vecHistEnergy.append(h_etaphiEnergy4)
vecHistEnergy.append(h_etaphiEnergy5)
vecHistEnergy.append(h_etaphiEnergy6)
vecHistEnergy.append(h_etaphiEnergy7)
vecHistEnergy.append(h_etaphiEnergy8)
vecHistEnergy.append(h_etaphiEnergy9)
vecHistEnergy.append(h_etaphiEnergy10)
vecHistEnergy.append(h_etaphiEnergy11)
vecHistEnergy.append(h_etaphiEnergy12)
vecHistEnergy.append(h_etaphiEnergy13)
vecHistEnergy.append(h_etaphiEnergy14)
vecHistEnergy.append(h_etaphiEnergy15)
vecHistEnergy.append(h_etaphiEnergy16)
vecHistEnergy.append(h_etaphiEnergy17)
vecHistEnergy.append(h_etaphiEnergy18)

minE = 1e-4
titles = []
titlesH = []
legendsE = []
legendsH = []
minX = -.3
maxX = .3
minY = -.3
maxY = .3

###################################################DRAW ENERGY MAP
maxEecal = 10. #h_etaphiEnergy15.GetMaximum()
maxEhcal = 10. #h_etaphiEnergy15.GetMaximum()

c1Energy = r.TCanvas("c1Energy","c1Energy",600,900)
r.gStyle.SetTitleSize(0.1,"xyz")
r.gStyle.SetLabelSize(0.1,"xyz")
c1Energy.Divide(3,6)

maximum = vecHistEnergy[9].GetMaximum()
minimum = vecHistEnergy[0].GetMinimum()

for i in range(1,19):
    titles.append('EMB layer '+str(i))
    titlesH.append('HB layer '+str(i))
    if i < 9:
        legendsE.append(r.TLegend(0.25,0.75,0.85,0.98, titles[i-1]))
        legendsE[i-1].SetTextSize(0.13)
        legendsE[i-1].SetFillColor(0)
        legendsE[i-1].SetFillStyle(0)
    else :
        legendsH.append(r.TLegend(0.25,0.75,0.85,0.98, titlesH[i-9]))
        legendsH[i-9].SetTextSize(0.13)
        legendsH[i-9].SetFillColor(0)
        legendsH[i-9].SetFillStyle(0)

    c1Energy.cd(i)
    c1Energy.cd(i).SetRightMargin(0.18)
    c1Energy.cd(i).SetLogz()
    if i==16:
        draw_hist2d(vecHistEnergy[i-1], "#bf{#Delta#eta}", "#bf{#Delta#phi}", "#bf{hit cell}", titles[i-1])
    else:
        draw_hist2d(vecHistEnergy[i-1], "", "", "", titles[i-1])
    
    vecHistEnergy[i-1].GetXaxis().SetNdivisions(4)
    vecHistEnergy[i-1].GetYaxis().SetNdivisions(4)
    vecHistEnergy[i-1].GetZaxis().SetNdivisions(3)
    vecHistEnergy[i-1].GetXaxis().SetRangeUser(minX,maxX)
    vecHistEnergy[i-1].GetYaxis().SetRangeUser(minY,maxY)
    vecHistEnergy[i-1].GetXaxis().SetTitleSize(0.15)
    vecHistEnergy[i-1].GetYaxis().SetTitleSize(0.15)
    vecHistEnergy[i-1].GetZaxis().SetTitleSize(0.15)
    vecHistEnergy[i-1].GetXaxis().SetLabelSize(0.12)
    vecHistEnergy[i-1].GetYaxis().SetLabelSize(0.12)
    vecHistEnergy[i-1].GetZaxis().SetLabelSize(0.15)
    vecHistEnergy[i-1].GetYaxis().SetTitleOffset(.7)
    vecHistEnergy[i-1].GetXaxis().SetTitleOffset(.7)
    vecHistEnergy[i-1].GetZaxis().SetTitleOffset(.5)
    vecHistEnergy[i-1].SetMaximum(maximum)         
    vecHistEnergy[i-1].SetMinimum(minimum)
    if i==16:
        vecHistEnergy[i-1].GetZaxis().SetLabelSize(0.1)

    if i < 9:
        legendsE[i-1].Draw()
    else:
        legendsH[i-9].Draw()

    c1Energy.cd(i).Modified()
    r.gPad.Modified()
    
    print 'Set Maximum: ', maximum

c1Energy.Update()
c1Energy.Write()
c1Energy.Print("rootFiles/bFieldOn/official/cells_resegmentedHCal_pi-"+str(int(Energy))+"GeV_energyMap.png")

#for i in range(1,19):
#    if i < 9:
#        th2f_to_csv(vecHistEnergy[i-1], "excelFiles/cells_scaled_log_layer"+str(i-1)+".csv", 1, etaBinsECal, 1, phiBinsECal, maximum, True)
#
#    else:
#        th2f_to_csv(vecHistEnergy[i-1], "excelFiles/cells_scaled_log_layer"+str(i-1)+".csv", 1, etaBinsHCal, 1, phiBinsHCal, maximum, True)


h_etaphiEnergy1.Write()
h_etaphiEnergy2.Write()
h_etaphiEnergy3.Write()
h_etaphiEnergy4.Write()
h_etaphiEnergy5.Write()
h_etaphiEnergy6.Write()
h_etaphiEnergy7.Write()
h_etaphiEnergy8.Write()
h_etaphiEnergy9.Write()
h_etaphiEnergy10.Write()
h_etaphiEnergy12.Write()
h_etaphiEnergy12.Write()
h_etaphiEnergy13.Write()
h_etaphiEnergy14.Write()
h_etaphiEnergy15.Write()
h_etaphiEnergy16.Write()
h_etaphiEnergy17.Write()
h_etaphiEnergy18.Write()

f.Close()

raw_input("Press ENTER to exit")
