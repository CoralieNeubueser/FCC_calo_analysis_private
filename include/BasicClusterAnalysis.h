#ifndef COMBINEDBASICCLUSTERANALYSIS_H
#define COMBINEDBASICCLUSTERANALYSIS_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class BasicClusterAnalysis: public BaseAnalysis  {

 public:
  BasicClusterAnalysis(double aEnergy, const std::string& aBitfieldEcal, const std::string& aClusterCollection, const std::string& aEcalCollection);
  ~BasicClusterAnalysis();

  Decoder ecal_decoder;
  Decoder* hcal_decoder = new Decoder("system:4,module:7,row:9,layer:5");
  std::string ecal_collection;
  std::string clusterCollection;

  TH1F* h_clusters;
  TH1F* h_nClusters;
  TH1F* h_ET;
  TH1F* h_phiRec;
  TH1F* h_etaRec;
  TH1F* h_thetaRec;
  TH1F* h_clusterEnergy;
  TH1F* h_clusterPt;
  TH1F* h_clusterCells;
  TH1F* h_cellId;

  TH2F* h_corrEnergyNumber;

  TH1F* h_ene_x;
  TH1F* h_ene_y;
  TH1F* h_ene_z;

  TH1F* h_ene_eta;
  TH1F* h_ene_phi;
  TH1F* h_ene_r;
  TH2F *h_lostECorr;
  TH1F* h_lambda;
  

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  double m_energy;
  double SumECluster;   // Cluster energy
  double SumET;    // Total transverse energy per event
 

};

#endif /* COMBINEDBASICCLUSTERANALYSIS_H */
