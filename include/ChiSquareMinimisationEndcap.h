#ifndef CHISQUAREMINIMISATIONENDCAP_H
#define CHISQUAREMINIMISATIONENDCAP_H

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

class ChiSquareMinimisationEndcap: public MultiFileAnalysis {

 public:
  ChiSquareMinimisationEndcap(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal);
  ~ChiSquareMinimisationEndcap();

  virtual void processEvent(podio::EventStore& store, int aEventId, int aEnergy, bool verbose) final;
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
  
  // Parameters from minimisation 5TeV 
  double a = 1.023;
  double b = 1.;
  double c = 0.446;
  double d = -1.73e-5;

  // Sampling fracions for EM scale
  double sf_ecal = .188; 
  double pion_hcal = .024;
  double sf_hcal = .027;

  double m_thr = 1e-6;
  double E_firstLayer;
  double E_lastLayer;
  double m_energy;

  double SumE_cell;    // Total cell energy per event in ECAL
  double SumH_cell;    // Total cell energy per event in HCAL
  
};

#endif /* CHISQUAREMINIMISATIONBARREL_H */
