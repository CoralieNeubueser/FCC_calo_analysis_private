#ifndef CHISQUAREMINIMISATIONPbSpacer_H
#define CHISQUAREMINIMISATIONPbSpacer_H

#include "BaseAnalysis.h"

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

class ChiSquareMinimisationPbSpacer: public BaseAnalysis {

 public:
  ChiSquareMinimisationPbSpacer(double aEnergy, std::vector<int> aE, std::vector<double> aVecEem, std::vector<double> aVecEhad, std::vector<double> aVecEfirst, std::vector<double> aVecElast);
  ~ChiSquareMinimisationPbSpacer();

  virtual void processEvent(podio::EventStore& store, int aEventId, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
 
 private:
  void Initialize_histos();
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
 
  // Parameters from minimisation 5TeV
  double a = 0.760792;//1.0716; //3.39;
  double b = 1.;//0.638; //3.39;
  double c = 0.978363;//-1.74e-05;//-6.4e-05;
  double d = 0.000437673;

  // Sampling fracions for EM scale
  double sf_ecal = .188; 
  double pion_hcal = .024;
  double sf_hcal = .027;

  double m_thr = 1e-6;
  double E_firstLayer;
  double E_lastLayer;

  std::vector<int>  m_E;
  std::vector<double> m_vecEem;
  std::vector<double> m_vecEhad;
  std::vector<double> m_vecEfirst;
  std::vector<double> m_vecElast; 
  double m_energy;

  double SumE_cell;    // Total cell energy per event in ECAL
  double SumH_cell;    // Total cell energy per event in HCAL

  const double RcaloMax_ecal = 2570.; // max ECAL, Cryostat 10cm + 5cm LAr bath decreased radius from 2750 to 2600
  const double RcaloMin_hcal = 2860.; // min HCAL
  const double layerThickness = 100; // mainly used as first layer thickness, 6 ECAL
  
};

#endif /* CHISQUAREMINIMISATIONPbSpacer_H */
