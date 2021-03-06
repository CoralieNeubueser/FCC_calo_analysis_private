#ifndef COMBINEDCLUSTERANALYSIS_H
#define COMBINEDCLUSTERANALYSIS_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ClusterAnalysis: public BaseAnalysis  {

 public:
  ClusterAnalysis(double aEnergy, const std::string& aBitfieldEcal, const std::string& aBitfieldHcal, const std::string& aClusterCollection, const std::string& aEcalCollection, bool normalizeToMCParticle, double deltaRCut);
  ~ClusterAnalysis();

  bool m_norm;
  double m_deltaRCut;
  
  Decoder ecal_decoder;
  Decoder hcal_decoder; 
  std::string ecal_collection;
  std::string clusterCollection;

  TH1F* h_clusters;
  TH1F* h_nClusters;
  TH1F* h_sharedClusters;
  TH1F* h_numCells;
  TH1F* h_ET;
  TH1F* h_phiRec;
  TH1F* h_etaRec;
  TH1F* h_thetaRec;
  TH1F* h_clusterEnergy;
  TH1F* h_clusterCells;
  TH1F* h_clusterCellTypes;
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
  TH2F *h_etaphi13;
  TH2F *h_etaphi14;
  TH2F *h_etaphi15;
  TH2F *h_etaphi16;
  TH2F *h_etaphi17;
  TH2F *h_etaphi18;

  /// eta from hits' direction
  TH2F *h_etaphiSharedCluster1;
  TH2F *h_etaphiSharedCluster2;
  TH2F *h_etaphiSharedCluster3;
  TH2F *h_etaphiSharedCluster4;
  TH2F *h_etaphiSharedCluster5;
  TH2F *h_etaphiSharedCluster6;
  TH2F *h_etaphiSharedCluster7;
  TH2F *h_etaphiSharedCluster8;
  TH2F *h_etaphiSharedCluster9;
  TH2F *h_etaphiSharedCluster10;
  TH2F *h_etaphiSharedCluster11;
  TH2F *h_etaphiSharedCluster12;
  TH2F *h_etaphiSharedCluster13;
  TH2F *h_etaphiSharedCluster14;
  TH2F *h_etaphiSharedCluster15;
  TH2F *h_etaphiSharedCluster16;
  TH2F *h_etaphiSharedCluster17;
  TH2F *h_etaphiSharedCluster18;

  /// eta from hits' direction                                                                                       
  TH2F *h_etaphiEnergy1;
  TH2F *h_etaphiEnergy2;
  TH2F *h_etaphiEnergy3;
  TH2F *h_etaphiEnergy4;
  TH2F *h_etaphiEnergy5;
  TH2F *h_etaphiEnergy6;
  TH2F *h_etaphiEnergy7;
  TH2F *h_etaphiEnergy8;
  TH2F *h_etaphiEnergy9;
  TH2F *h_etaphiEnergy10;
  TH2F *h_etaphiEnergy11;
  TH2F *h_etaphiEnergy12;
  TH2F *h_etaphiEnergy13;
  TH2F *h_etaphiEnergy14;
  TH2F *h_etaphiEnergy15;
  TH2F *h_etaphiEnergy16;
  TH2F *h_etaphiEnergy17;
  TH2F *h_etaphiEnergy18;
  /// eta from hits' direction                                                                                       
  TH2F *h_etaphiType1;
  TH2F *h_etaphiType2;
  TH2F *h_etaphiType3;
  TH2F *h_etaphiType4;
  TH2F *h_etaphiType5;
  TH2F *h_etaphiType6;
  TH2F *h_etaphiType7;
  TH2F *h_etaphiType8;
  TH2F *h_etaphiType9;
  TH2F *h_etaphiType10;
  TH2F *h_etaphiType11;
  TH2F *h_etaphiType12;
  TH2F *h_etaphiType13;
  TH2F *h_etaphiType14;
  TH2F *h_etaphiType15;
  TH2F *h_etaphiType16;
  TH2F *h_etaphiType17;
  TH2F *h_etaphiType18;

  THStack *s_EPhi;

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

  double SumECluster;   // Cluster energy
  double SumET;    // Total transverse energy per event
 
  // hcal sampling fraction correction to .0774%
  double hcalCorr_sf = 0.9299;//0.00077/0.000828;
  double lambdaOffset = 0.2; //#lambda
  double lambdaOffsetHcal = 3.; //#lambda
  double zMinEcal = 16640; //*mm
  double lambdaEcal = 164.2; //mm
  double lambdaHcal = 156.1; //mm
  double zMinHcal = 17120; //*mm

};

#endif /* COMBINEDCLUSTERANALYSIS_H */
