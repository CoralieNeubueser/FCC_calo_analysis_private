#ifndef CELLSIGNIFICANCES_H
#define CELLSIGNIFICANCES_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CellSignificances: public BaseAnalysis  {

 public:
  CellSignificances(double aEnergy, const std::string& aBitfieldEcal, const std::string& aHitCollection, bool aElectrNoiseAdded);
  ~CellSignificances();

  bool elecNoise;
  Decoder ecal_decoder;
  Decoder* hcal_decoder = new Decoder("system:4,module:7,row:9,layer:5");
  std::string hit_collection;
  
  TH1F* h_energy_EMlayer1;
  TH1F* h_energy_EMlayer2;
  TH1F* h_energy_EMlayer3;
  TH1F* h_energy_EMlayer8;
  
  TH1F* h_energy_HadLayer1;
  TH1F* h_energy_HadLayer2;
  TH1F* h_energy_HadLayer6;
  TH1F* h_energy_HadLayer10;
  
 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  void Initialize_histos();

  double m_energy;
  double m_thr_ecal = 0.0075/4.;
  // noise per layer at ca. 0.36
  double m_thr_ecal_layer1 = 0.004;
  double m_thr_ecal_layer2 = 0.01;
  double m_thr_ecal_layer3 = 0.0045;
  double m_thr_ecal_layer8 = 0.015;

  double m_thr_hcal = 0.0115/4.;
  double m_thr_hcal_tile = 0.009;

  std::map<uint64_t, std::pair<double,double>> m_map;
};

#endif /* CELLSIGNIFICANCES_H */
