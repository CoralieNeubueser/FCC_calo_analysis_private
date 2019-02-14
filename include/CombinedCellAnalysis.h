#ifndef COMBINEDCELLANALYSIS_H
#define COMBINEDCELLANALYSIS_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedCellAnalysis: public BaseAnalysis  {

 public:
  CombinedCellAnalysis(double aEnergy, double a1, double a2, double a3, double b, double c1, double c2, double c3, double d, double aLinCorrA, double aLinCorrB, double aLinCorrC, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal);
  ~CombinedCellAnalysis();

  Decoder ecal_decoder;
  Decoder hcal_decoder;

  TH1F* h_benchmark;
  TH1F* h_benchmarkLinCorr;
  TH1F* h_cellEnergy;
  TH1F* h_numCells;
  TH1F* h_cellEnergy_ecal;
  TH1F* h_cellEnergy_hcal;

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  // Sampling fracions for EM scale
  double sf_ecal = .188; 
  double pion_hcal = .029;
  double sf_hcal = .032;

  double m_energy;
  double m_thr;
  double m_a1;
  double m_a2;
  double m_a3;
  double m_b; 
  double m_c1;
  double m_c2;
  double m_c3;
  double m_d;

  double m_linCorrA = 0.;
  double m_linCorrB = 1.;
  double m_linCorrC = 1.;
 
  double SumE_ecal;    // Total hit energy per event
  double SumE_hcal;    // Total hit energy per event
  double E_firstLayer;    // Energy in first HCAL layer calibrated to pions
  double E_lastLayer;    // Energy in last ECAl layer on EM scale 

  const double RcaloMax_ecal = 2600.; // max ECAL, Cryostat 10cm + 5cm LAr bath decreased radius from 2750 to 2600
  const double RcaloMin_hcal = 2860.; // min HCAL
  const double layerThickness = 100; // mainly used as first layer thickness, 6 ECAL

  // Parameters from minimisation 5TeV 
  // double a = 1.0716; //3.39;
  // double b = 0.638; //3.39;
  // double c = -1.74e-05;//-6.4e-05;                                            

};

#endif /* COMBINEDCELLANALYSIS_H */
