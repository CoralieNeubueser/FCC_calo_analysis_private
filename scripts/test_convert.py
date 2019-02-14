from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import CombinedCellAnalysis, ClusterAnalysis, CombinedCellAnalysisPbSpacer, TCanvas, TFile, TF1, TH2D, THStack, gPad, TH1D, TMath, TPad, gROOT, gStyle, TStyle, TLegend
from draw_functions import draw_1histogram, draw_2histograms, draw_hist2d
from export_to_csv import *

gROOT.SetBatch(True)

Energy = 100

import glob
# /eos/experiment/fcc/hh/simulation/samples/v03/singlePart/pim/bFieldOff/eta0.36/10GeV/ntup/positions/resegmentedHCal/electronicsNoise/coneCut/0.4/noSignal/output_22113320.root
File = '/eos/experiment/fcc/hh/simulation/samples/v03/singlePart/pim/bFieldOff/eta0.36/'+str(int(Energy))+'GeV/ntup/positions/resegmentedHCal/electronicsNoise/coneCut/0.4/noSignal/output_*.root'
print File 

for i,ifile in enumerate(glob.glob(File)):
    print i
