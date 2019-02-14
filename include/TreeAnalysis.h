#ifndef TREEANALYSIS_H
#define TREEANALYSIS_H

#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include <string>

namespace podio {
  class EventStore;
  class ROOTReader;
}

class TreeAnalysis {
public:
  TreeAnalysis();
  virtual ~TreeAnalysis();

  void loop(const std::string& aFilename, bool aVerbose = false);  //Open the file in the reader and loop through the events
  protected:
  TTree* m_tree;
  TFile* m_file;

private:
  virtual void processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose = false) = 0;
  virtual void finishLoop(int aNumEvents, bool aVerbose) = 0;
};

#endif /* TREEANALYSIS_H */
