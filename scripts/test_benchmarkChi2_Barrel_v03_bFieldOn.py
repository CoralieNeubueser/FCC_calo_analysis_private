# configure (get file name etc.)
import calo_init
calo_init.add_defaults()
calo_init.parser.add_argument("--bitfieldEcal", help="bitfield of dd4hep ECAL readout", type = str)
calo_init.parser.add_argument("--bitfieldHcal", help="bitfield of dd4hep HCAL readout", type = str)
calo_init.parse_args()
calo_init.print_config()

#Setup ROOT
from ROOT import gSystem
gSystem.Load("libCaloAnalysis")
from ROOT import ChiSquareMinimisationBarrel, gStyle, TCanvas, TFile, TH1F, TH2F, TF1, TProfile, gPad, TMath, THStack, TGraphErrors
import numpy as np
from array import array
import ctypes
#import draw functions
from draw_functions import draw_1histogram, draw_2histograms

Energies = [10,100,1000,10000]
ecalBitfield = "system:4,cryo:1,type:3,subtype:3,layer:8,eta:9,phi:10"
hcalBitfield = "system:4,module:8,row:9,layer:5"

if calo_init.args.bitfieldEcal:
    ecalBitfield = calo_init.args.bitfieldEcal
if calo_init.args.bitfieldHcal:
    hcalBitfield = calo_init.args.bitfieldHcal

h_parameter = TH1F("h_parameter","", 5, 0,5);
h_energyECAL = TH1F("h_energyECAL","", 100, 0, Energies[0]);
h_energyHCAL = TH1F("h_energyHCAL","", 100, 0, Energies[0]);
h_benchmark = TH1F("h_benchmark","", 1000, 0, Energies[0]+Energies[0]*0.5);
h_emScale = TH1F("h_emScale","", 1000, 0, Energies[0]+Energies[0]*0.5);

meanA =0.
meanB =0.
meanC =0.
numFiles=0

files = []
energies = []

import glob
for iEnergy in Energies:
    i = 0
    for ifile in glob.glob('/eos/experiment/fcc/hh/simulation/samples/v03/singlePart/pim/bFieldOn/eta0.36/'+str(iEnergy)+'GeV/simu/*.root'): 
        if i > 5: 
            break
        print ifile
        files.append(ifile)
        energies.append(iEnergy)
        i += 1
print len(files)
vecF = np.array(files)
arrF = (ctypes.c_char_p * len(files))(*files)
arrE = (ctypes.c_double * len(energies))(*energies)
numF = len(arrF)

ma = ChiSquareMinimisationBarrel(float(Energies[0]), ecalBitfield, hcalBitfield)
ma.loop(int(numF), arrF, arrE, calo_init.verbose)
print "ChiSquare Parameters: a = ", ma.h_parameter.GetBinContent(0) ," +/- ", ma.h_parameter.GetBinError(0)  
print "ChiSquare Parameters: b = ", ma.h_parameter.GetBinContent(1) ," +/- ", ma.h_parameter.GetBinError(1)  
print "ChiSquare Parameters: c = ", ma.h_parameter.GetBinContent(2) ," +/- ", ma.h_parameter.GetBinError(2)  
print "ChiSquare Parameters: d = ", ma.h_parameter.GetBinContent(3) ," +/- ", ma.h_parameter.GetBinError(3)  

h_parameter.Add(ma.h_parameter,1)
h_energyECAL.Add(ma.h_energyECAL,1)
h_energyHCAL.Add(ma.h_energyHCAL,1)
h_benchmark.Add(ma.h_benchmark,1)
h_emScale.Add(ma.h_emScale,1)

print "Resolution EM scale:  ", h_emScale.GetFunction("gaus").GetParameter(2)/h_emScale.GetFunction("gaus").GetParameter(1)*100  
print "Resolution benchmark: ", h_benchmark.GetFunction("gaus").GetParameter(2)/h_benchmark.GetFunction("gaus").GetParameter(1)*100  

raw_input("Press ENTER to exit")
#c1.SaveAs("plots_"+PARTICLE+str(ENERGY)+".gif")
