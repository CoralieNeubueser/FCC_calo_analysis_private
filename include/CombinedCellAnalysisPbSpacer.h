#ifndef COMBINEDCELLANALYSISPBSPACER_H
#define COMBINEDCELLANALYSISPBSPACER_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CombinedCellAnalysisPbSpacer: public BaseAnalysis  {

 public:
  CombinedCellAnalysisPbSpacer(double aEnergy, double aThr, double a, double b, double c);
  CombinedCellAnalysisPbSpacer(double aEnergy, double aThr, double a, double b, double c, double aLin, double bLin, double cLin, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal);
  ~CombinedCellAnalysisPbSpacer();

  Decoder ecal_decoder;
  Decoder hcal_decoder;

  TH1D* h_benchmark;
  TH1D* h_benchmark2nd;
  TH1D* h_cellEnergy;
  TH1D* h_cellEnergy_ecal;
  TH1D* h_cellEnergy_hcal;
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
  TH2D *h_lostECorr;
  TH2D *h_lostEMultiCorr;
  
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

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  double m_energy;
  double m_thr;
  double m_a;
  double m_b; 
  double m_c;
  double m_aLin;
  double m_bLin; 
  double m_cLin;
 
  double SumE_ecal;    // Total hit energy per event
  double SumE_hcal;    // Total hit energy per event
  double E_firstLayer;    // Energy in first HCAL layer calibrated to pions
  double E_lastLayer;    // Energy in last ECAl layer on EM scale 

  // hcal sampling fraction correction to 2.98%
  double sf_hcalCorr = 1.;//(1./34.72)/0.0298;
  double lambdaOffset = 0.3; //#lambda
  double lambdaOffsetHcal = 2.2; //#lambda
  double rMinEcal = 1920; //*mm
  double lambdaEcal = 294; //mm
  double lambdaHcal = 200.7; //mm
  double rMinHcal = 2860.5; //*mm
};

#endif /* COMBINEDCELLANALYSISPBSPACER_H */
