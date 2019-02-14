#ifndef COMBINEDCELLANALYSISFWD_H
#define COMBINEDCELLANALYSISFWD_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedCellAnalysisFwd: public BaseAnalysis  {

 public:
  CombinedCellAnalysisFwd(double aEnergy, double aThr, double a, double b, double c, double d);
  CombinedCellAnalysisFwd(double aEnergy, double aThr, double a, double b, double c, double d, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aEcalCollection, const std::string& aHcalCollection);
  ~CombinedCellAnalysisFwd();

  Decoder ecal_decoder;
  Decoder hcal_decoder;
  std::string ecal_collection;
  std::string hcal_collection;

  TH1F* h_benchmark;
  TH1F* h_cellEnergy;
  TH1F* h_cellEnergy_ecal;
  TH1F* h_cellEnergy_hcal;
  TH1F* h_sf_ecal;
  TH1F* h_sf_hcal;
  TH1F* h_ET;
  TH1F* h_ET_ecal;
  TH1F* h_ET_hcal;
  TH1F* h_phiRec_ecal;
  TH1F* h_phiRec_hcal;
  TH1F* h_phiRec;
  TH1F* h_etaRec;
  TH1F* h_thetaRec;
  TH1F* h_etaRec_ecal;
  TH1F* h_etaRec_hcal;
  TH1F* h_thetaRec_ecal;
  TH1F* h_thetaRec_hcal;
  TH1F* h_cellEnergy_ecalP;
  TH1F* h_cellEnergy_hcalP;
  TH1F* h_cellEnergy_first;
  TH1F* h_cellEnergy_last;
  TH1F* h_cellId;

  TH1F* h_ene_x;
  TH1F* h_ene_y;
  TH1F* h_ene_z;

  TH1F* h_ene_eta;
  TH1F* h_ene_phi;
  TH1F* h_ene_r;
  TH2F *h_lostECorr;
  TH1F* h_lambda;
  
  /// eta from hits' direction
  TH2F *h_etaphi1;
  TH2F *h_etaphi2;
  TH2F *h_etaphi3;
  TH2F *h_etaphi4;
  TH2F *h_etaphi5;
  TH2F *h_etaphi6;
  TH2F *h_etaphi7;
  TH2F *h_etaphi8;
  TH2F *h_etaphi9;
  TH2F *h_etaphi10;
  TH2F *h_etaphi11;
  TH2F *h_etaphi12;

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
  double SumET_ecal;    // Total transverse energy per event
  double SumET_hcal;    // Total transverse energy per event
  double E_firstLayer;    // Energy in first HCAL layer calibrated to pions
  double E_lastLayer;    // Energy in last ECAl layer on EM scale 

  // hcal sampling fraction correction to .0774%
  double hcalCorr_sf = 0.9299;//0.00077/0.000828;
  double lambdaOffset = 0.2; //#lambda
  double lambdaOffsetHcal = 3.; //#lambda
  double zMinEcal = 16640; //*mm
  double lambdaEcal = 164.2; //mm
  double lambdaHcal = 156.1; //mm
  double zMinHcal = 17120; //*mm

};

#endif /* COMBINEDCELLANALYSISPBSPACER_H */
