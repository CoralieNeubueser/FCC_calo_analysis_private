#ifndef COMBINEDCELLANALYSISENDCAP_H
#define COMBINEDCELLANALYSISENDCAP_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedCellAnalysisEndcap: public BaseAnalysis  {

 public:
  CombinedCellAnalysisEndcap(double aEnergy, double aThr, double a, double b, double c, double d);
  CombinedCellAnalysisEndcap(double aEnergy, double aThr, double a, double b, double c, double d, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aBitfieldHcalEB, const std::string& aEcalCollection, const std::string& aHcalCollection, const std::string& aHcalEBCollection);
  ~CombinedCellAnalysisEndcap();

  Decoder ecal_decoder;
  Decoder hcal_decoder;
  Decoder hcalEB_decoder;
  std::string ecal_collection;
  std::string hcal_collection;
  std::string hcalEB_collection;

  TH1D* h_benchmark;
  TH1D* h_cellEnergy;
  TH1D* h_ET;
  TH1D* h_cellEnergy_ecal;
  TH1D* h_cellEnergy_hcal;
  TH1D* h_cellEnergy_hcalEB;
  TH1D* h_sf_ecal;
  TH1D* h_sf_hcal;
  TH1D* h_ET_ecal;
  TH1D* h_ET_hcal;
  TH1D* h_ET_hcalEB;
  TH1D* h_phiRec_ecal;
  TH1D* h_phiRec_hcal;
  TH1D* h_phiRec;
  TH1D* h_etaRec;
  TH1D* h_thetaRec;
  TH1D* h_etaRec_ecal;
  TH1D* h_etaRec_hcal;
  TH1D* h_thetaRec_ecal;
  TH1D* h_thetaRec_hcal;
  TH1D* h_cellEnergy_ecalP;
  TH1D* h_cellEnergy_hcalP;
  TH1D* h_cellEnergy_first;
  TH1D* h_cellEnergy_last;
  TH1D* h_cellId;

  TH1D* h_ene_x;
  TH1D* h_ene_y;
  TH1D* h_ene_z;

  TH1D* h_ene_eta;
  TH1D* h_ene_phi;
  TH1D* h_ene_r;
  TH1D* h_lambdaEcal;
  TH1D* h_lambdaHcal;
  TH1D* h_lambdaHcalEB;
  TH2D *h_lostECorr;
  
  /// eta from hits' direction
  TH2D *h_etaphi1;
  TH2D *h_etaphi2;
  TH2D *h_etaphi3;
  TH2D *h_etaphi4;
  TH2D *h_etaphi5;
  TH2D *h_etaphi6;

  TH2D *h_etaphi7;
  TH2D *h_etaphi8;
  TH2D *h_etaphi9;
  TH2D *h_etaphi10;
  TH2D *h_etaphi11;
  TH2D *h_etaphi12;

  TH2D *h_etaphi13;
  TH2D *h_etaphi14;
  TH2D *h_etaphi15;
  TH2D *h_etaphi16;
  TH2D *h_etaphi17;
  TH2D *h_etaphi18;
  TH2D *h_etaphi19;
  TH2D *h_etaphi20;

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  double m_energy;
  double m_thr;
  double m_a;
  double m_b; 
  double m_c;
  double m_d;

  double SumE_ecal;    // Total hit energy per event
  double SumE_hcal;    // Total hit energy per event
  double SumE_hcalEB;    // Total hit energy per event
  double SumET_ecal;    // Total transverse energy per event
  double SumET_hcal;    // Total transverse energy per event
  double SumET_hcalEB;    // Total transverse energy per event
  double E_firstLayer;    // Energy in first HCAL layer calibrated to pions
  double E_lastLayer;    // Energy in last ECAl layer on EM scale 

  // hcal sampling fraction correction to 2.98%
  double hcalCorr_sf = 1.0049;//0.02988/0.02974;
  double lambdaOffset = 0.2; //#lambda
  double lambdaOffsetHcal = 3.; //#lambda
  double zMinEcal = 5440; //*mm
  double lambdaEcal = 277.2; //mm
  double lambdaHcal = 178.7; //mm
  double zMinHcal = 6020; //*mm
};

#endif /* COMBINEDCELLANALYSISENDCAP_H */
