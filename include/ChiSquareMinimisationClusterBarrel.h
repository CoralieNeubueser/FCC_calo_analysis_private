#ifndef CHISQUAREMINIMISATIONCLUSTERBARREL_H
#define CHISQUAREMINIMISATIONCLUSTERBARREL_H

#include "MultiFileAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "Math/IFunction.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ChiSquareMinimisationClusterBarrel: public MultiFileAnalysis {

 public:
  ChiSquareMinimisationClusterBarrel(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, bool aEnergyDependence);
  ~ChiSquareMinimisationClusterBarrel();

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

  /// Energy measured in ECAL on EM scale
  TH1F* h_energyEMcluster;
  /// Energy measured in HCAL non-calibrated 
  TH1F* h_energyHADcluster;

  /// Energy reconstruction on EM scale
  TH1F* h_emScale;
  /// Benchmark reconstructed energy
  TH1F* h_benchmark;
  
  Decoder* m_decoder = new Decoder("system:4");
  Decoder m_decoderECal;
  Decoder m_decoderHCal;
  bool m_eDep = false;
  
 private:
  void Initialize_histos();
  
  // Parameters from minimisation 100GeV & 5TeV 
  double a = 1.00537;//1.00537,1,0.406406,-1.39058e-050.96475,1,0.455803,-2.2142e-06
  double b = 1.;
  double c = 0.406406;//0.455803;
  double d = -1.39058e-05;//-2.2142e-06;

  double m_energy;
  double m_fractionECal=0.9;

  /// e/h of HCal
  double m_ehHCal = 1.1;

  int m_lastECalLayer = 7;
  int m_firstHCalLayer = 0;

  uint m_systemIdHCal = 8;
  uint m_systemIdECal = 5;
  
};

#endif /* CHISQUAREMINIMISATIONCLUSTERBARREL_H */
