#ifndef COMBINEDCELLANALYSISBENCHMARKBARREL_H
#define COMBINEDCELLANALYSISBENCHMARKBARREL_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedCellAnalysisBenchmarkBarrel: public BaseAnalysis  {

 public:
  CombinedCellAnalysisBenchmarkBarrel(double aEnergy, double aThr, double a, double b, double c, double aLin, double bLin, double cLin, double a0, double a1, double a2, double b0, double b1, double b2, double c0, double c1, double c2, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal);
  ~CombinedCellAnalysisBenchmarkBarrel();

  Decoder ecal_decoder;
  Decoder hcal_decoder;

  TH1D* h_classicBenchmark;
  TH1D* h_eDepBenchmark;
  TH1D* h_eTot;

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  double m_energy;
  double m_thr;
  double m_a;
  double m_b; 
  double m_c;
  double m_a0;
  double m_a1;
  double m_a2;
  double m_b0;
  double m_b1;
  double m_b2;
  double m_c0;
  double m_c1;
  double m_c2;
  double m_aLin;
  double m_bLin; 
  double m_cLin;
 
  double SumE_ecal;    // Total hit energy per event
  double SumE_hcal;    // Total hit energy per event
  double E_firstLayer;    // Energy in first HCAL layer calibrated to pions
  double E_lastLayer;    // Energy in last ECAl layer on EM scale 

};

#endif /* COMBINEDCELLANALYSISBENCHMARKBARREL_H */
