#ifndef CHISQUAREMINIMISATIONEDEPBARREL_H
#define CHISQUAREMINIMISATIONEDEPBARREL_H

#include "MultiFileAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "Math/IFunction.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ChiSquareMinimisationEdepBarrel: public MultiFileAnalysis {

 public:
  ChiSquareMinimisationEdepBarrel(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal);
  ~ChiSquareMinimisationEdepBarrel();

  virtual void processEvent(podio::EventStore& store, int aEventId, double aEnergy, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;

  /// Results of minimisation 
  // bin=0: HCAL calibration factor for hadrons
  // bin=1: correction factor for cryostat
  // bin=2: corection factor for non-compensating ECAL 
  TH1F* h_parameter;

  /// Energy measured in ECAL on EM scale
  TH1F* h_energyECAL;
  /// Energy measured in HCAL non-calibrated
  TH1F* h_energyHCAL;
  /// Energy reconstruction on EM scale
  TH1F* h_emScale;
  /// Benchmark reconstructed energy
  TH1F* h_benchmark;
  
  Decoder ecal_decoder;
  Decoder hcal_decoder;

 private:
  void Initialize_histos();
  
  // Parameters from minimisation 100GeV & 5TeV 
  double a = 1.00537;//1.00537,1,0.406406,-1.39058e-050.96475,1,0.455803,-2.2142e-06
  double b = 1.;
  double c = 0.406406;//0.455803;
  double d = -1.39058e-05;//-2.2142e-06;

  // Sampling fracions for EM scale
  double sf_ecal = .188; 
  double pion_hcal = .024;
  double sf_hcal = .027;

  double m_thr = 1e-6;
  double EHCAL_firstLayer;
  double EECAL_firstLayer;
  double EECAL_lastLayer;
  double m_energy;

  double SumE_cell;    // Total cell energy per event in ECAL
  double SumEcorr_cell;// Corrected for upstream material ECAL
  double SumH_cell;    // Total cell energy per event in HCAL
 
};

#endif /* CHISQUAREMINIMISATIONEDEPBARREL_H */
