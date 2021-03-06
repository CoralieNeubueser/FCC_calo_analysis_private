#ifndef COMBINEDCLUSTERANALYSISTREE_H
#define COMBINEDCLUSTERANALYSISTREE_H

#include "TreeAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ClusterAnalysisTree: public TreeAnalysis  {

 public:
  ClusterAnalysisTree(double aEnergy, const std::string& aBitfieldEcal, const std::string& aClusterCollection, const std::string& aEcalCollection, const std::string& aRootFileName);
  ~ClusterAnalysisTree();

 private:

  virtual void processEvent(podio::EventStore& store, int aEventId, bool aVerbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  
  void Initialize_tree();
  void clearVariables();

  TFile* m_file;

  Decoder ecal_decoder;
  std::string clusterCollection;
  std::string ecal_collection;
  std::string m_rootFileName;
  
  int m_event = 5;
  /// tree variables
  double m_Energy;
  int m_cluster;
  std::vector<double> m_clusterEnergy;
  std::vector<int> m_clusterCells;
  std::vector<std::vector<double> >  m_etaCell;
  std::vector<std::vector<double> >  m_thetaCell;
  std::vector<std::vector<double> >  m_phiCell;
  std::vector<std::vector<int> >  m_layerCell;
  std::vector<std::vector<int> >  m_typeCell;
  std::vector<std::vector<double> >  m_energyCell;

  double m_energy;

  double SumECluster;   // Cluster energy
  double SumET;    // Total transverse energy per event
 
};

#endif /* COMBINEDCLUSTERANALYSISTREE_H */
