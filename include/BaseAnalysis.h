#ifndef BASEANALYSIS_H
#define BASEANALYSIS_H

#include "TObject.h"
#include "TH1.h"
#include "THStack.h"
#include "TProfile.h"
#include <string>

namespace podio {
  class EventStore;
  class ROOTReader;
}

class BaseAnalysis {
public:
  BaseAnalysis();
  virtual ~BaseAnalysis();

  void loop(const std::string& aFilename, bool aVerbose = false, bool aOneEvent = false);  //Open the file in the reader and loop through the events
  protected:
  std::vector<TH1*> m_histograms;
  std::vector<THStack*> m_stacks;
  std::vector<TProfile*> m_profiles;
private:
  virtual void processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose = false) = 0;
  virtual void finishLoop(int aNumEvents, bool aVerbose) = 0;
  virtual void Delete_histos();
  virtual void Reset_histos();
};

#endif /* BASEANALYSIS_H */
